package Xcalibr;

use strict;
use warnings;

use 5.010;

#use Data::Dumper; # qw/dump/;
#use Sort::Naturally qw/nsort/;
use Carp;
use SequenceSet;

sub new {
	my $class = shift;
	my $settings = shift;

	my $self = {};
	bless $self, $class;

	#TODO check settings hash for required keys
	$self->{settings} = $settings;

	return $self;
}

sub parse_template {
	my $self = shift;
	my $desc = shift;

    #croak "Error in sequence template: '$desc'\n" if $desc !~ m/(([XYZNCGAT]{1})(n?(?!$)|\d*))+/;
	my $pm = qr/
		   ([XYZCGATN]{1}) #find a base or index once
		   (               #the n or number of repeats
		      n$ |         # an n is the end of the pattern OR
		      (\d*)        #the number the base should repeat
		   )
		/x;
	croak "Error in sequence template: '$desc'\n" if $desc !~ m/^($pm)+$/x;
	my $pos=0;
	while ( $desc =~ m/$pm/g ) {
		my $type = $1;
		#if length is n the this will result in length 0 which will result in a clip to end substr function
		my $length = $2 eq "n" ? 0 : length($2) > 0 ? int($2) : 1;
		$self->{template}->{sequence} .= "$type" x $length;
		$self->{template}->{sequence} .= $type . "n" if $length == 0;

		given ($type) {
			when (/[XYZ]/) {
				die "Only one region per idx allowed (e.g. don't write ZZZ but Z3)" if exists $self->{template}{$_};
                $self->{template}{$_} = {
                    length    => $length,
                    pos       => $pos,
                    substrfun => make_substring_function(
                        $pos, $length, $self->{settings}->{reltoseq}),
                };
				continue;
			}
			default { $pos += $length; }
		}
	}
	(my $re = $self->{template}->{sequence}) =~ s/[^CGAT]/./g;
	#trim "." from start end
	$re =~ s/^\.+//;
	$re =~ s/\.+$//;

	#if $re is empty there is no fixed part to look for. There should probably be a minimum length 
	#or disable reltoseq
	if($self->{settings}->{reltoseq} && length($re) < 10) {
		croak "Fixed part is too small to enable relative clipping";
	}

	my $reg = qr/$re/;
	$self->{template}->{regexp} = $reg;
	#find expected template start by applying regex to self
	$self->{template}->{sequence} =~ m/$reg/;
	$self->{template}->{startpos} = $-[0];
}

sub make_substring_function {
	my $start = shift;
	my $length = shift;
	my $relative = shift;

	if ($length == 0) {
		return sub {my $str = shift; my $pos = $start+(shift // 0); $pos = 0 if $pos < 0; substr($str, $pos)} if $relative;
		return sub {substr(shift, $start)}
	} else {
		return sub {my $str = shift; my $pos = $start+(shift // 0); $pos = 0 if $pos < 0; substr($str, $pos, $length)} if $relative;
		return sub {substr(shift, $start, $length)};
	}
}

sub add_references {
	my $self = shift;
	my $refs = shift;
	my $mismatches = shift;
	foreach my $ref (keys %$refs) {
		if (defined $refs->{$ref}) {
			$self->add_part_referencefile($ref, $refs->{$ref}, $mismatches->{$ref});
		} elsif (defined $mismatches->{$ref}) {
			croak "--mismatches$ref requires a --match$ref <file.fa>";
		}
	}
	#TODO check for reference sequence for splitby
}

sub add_part_referencefile {
	my $self = shift;
	my $part = shift;
	my $file = shift;
	my $mismatches = shift;

	croak "Reference file $file for $part cannot be read." unless -e $file;

	#check for part in template
	croak "This part is not assigned in the template" unless exists $self->{template}->{$part};

	#load the file as sequenceset
	my $s = SequenceSet->new;
	$s->read_from_fasta($file, $mismatches);
	#add to self
	$self->{refsets}->{$part} = $s;
}

sub count_hash {
	my $self = shift;
	my $data = shift;

	croak "Invalid hashref" unless ref($data) eq "HASH";

	#we can now create the final match functions using the settings and the
	#loaded reference file(s)
	for my $t (qw/splitby rows cols/) {
		my $nameintemplate = $self->{settings}->{split}->{$t} ;
		#Did the user put this clause in the template?
		if(exists $self->{template}->{$nameintemplate}) {
			#the substr function is already generated. The match function is dependent on
			#the printall flag

			#if we have a reference try tyo find the name in the set
			if (exists $self->{refsets}->{$nameintemplate}) {
				if($self->{settings}->{printall}->{$nameintemplate}) {
					$self->{matchfuncs}->{$t} = sub {
						my $str = &{ $self->{template}->{$nameintemplate}->{substrfun} }(shift, shift);
						$self->{refsets}->{$nameintemplate}->match($str) // $str;
					};
				} else {
					$self->{matchfuncs}->{$t} = sub {
						$self->{refsets}->{$nameintemplate}->match(
							&{ $self->{template}->{$nameintemplate}->{substrfun} }(
								shift, shift ) ) // "nohit_$t";
					};
				}
			} else {
				#no ref so return the substring
				$self->{matchfuncs}->{$t} = $self->{template}->{$nameintemplate}->{substrfun};
			}
		} else {
			#this placeholder is unused.
			$self->{matchfuncs}->{$t} = sub { "nohit_$t" };
		}
	}


	#TODO implement printzero (by prefilling result structure)?
	for my $t (qw/splitby rows cols/) {
		my $nameintemplate = $self->{settings}->{split}->{$t} ;
		next unless defined $self->{settings}->{printempty}->{$nameintemplate};
		carp "WARNING: Printall requested but no reference supplied" && next unless defined $self->{refsets}->{$nameintemplate};
		$self->{countresults}->{seen}->{$t}->{$_} = 0 foreach $self->{refsets}->{$nameintemplate}->all_names;
	}

	my $reqseq = $self->{settings}->{requireseq};
	my $reltoseq = $self->{settings}->{reltoseq};
	foreach my $seq (keys %$data) {
		$self->{stats}->{'Total sequences seen'} += $data->{$seq};
		$self->{stats}->{'Total unique sequences seen'}++;
		my $offset = 0;
		#if the constant sequence is required test for it and optionally change match offset
		if($reqseq) {
			if($seq =~ /$self->{template}->{regexp}/) {
				#compare match pos with expected start
				$offset = $-[0] - $self->{template}->{startpos};
				if($reltoseq) {
					if($offset != 0) {
						$self->{stats}->{'Moved offset based on template sequence'}++;
						#FIXME if offset is larger than pos in substr fun this function clips of the
						#end of the string which leads to unpredictable results
					}
				}
				#fail sequence if template found but offset !=0 and relative is off
				if ($offset != 0 && (not defined $reltoseq)) {
					$self->{stats}->{'Template sequence found at wrong location'}++;
					next;
				}
			} else {
				$self->{stats}->{'Template sequence not found'}++;
				next;
			}
		}
		#match all regions from template
		my $splitby = &{$self->{matchfuncs}->{splitby}}($seq, $offset);
		my $row = &{$self->{matchfuncs}->{rows}}($seq, $offset);
		my $col = &{$self->{matchfuncs}->{cols}}($seq, $offset);
		
		#store seen tags/barcode/index stuff for easy printing (this wastes memory)
		$self->{countresults}->{seen}->{splitby}->{$splitby} = 0;
		$self->{countresults}->{seen}->{rows}->{$row} = 0;
		$self->{countresults}->{seen}->{cols}->{$col} = 0;

		#and link table results in hash(ofhashofhash)
		$self->{countresults}->{$splitby}->{$row}->{$col} +=  $data->{$seq};
	}
}

sub write_tables {
	my $self = shift;
	my $output = shift;

	#use nsort for small key sets ?? Maybe Sort::Key?
	my @split = sort keys %{$self->{countresults}->{seen}->{splitby}};
	my @rows  = sort keys %{$self->{countresults}->{seen}->{rows}};
	my @cols  = sort keys %{$self->{countresults}->{seen}->{cols}};


	#write to file or stdout (if not splitby and no --output)
	foreach my $splitfile (@split) {
		#FIXME  if somebody names it splitfile splitby we're in deep trouble....
		if ($output) {
			my $filename = $output;
			$filename =~ s/(\.txt)?$/-$splitfile.txt/i 
			   unless $splitfile eq "nohit_splitby" && not exists $self->{template}->{$self->{settings}->{split}->{splitby}};
			#add .txt to filename if necessary
			$filename .= "\.txt" unless $filename =~ /\.txt$/i;
			open(OUT,">",$filename) or croak "$filename:$!";
			select(OUT);
		}
		#print header
		say join("\t", "name", @cols);
		foreach my $row (@rows) {
			say join("\t", $row, map { $self->{countresults}->{$splitfile}->{$row}->{$_} // 0} @cols)
				if defined $self->{settings}->{printempty}->{$self->{settings}->{split}->{rows}} or
					scalar(keys %{$self->{countresults}->{$splitfile}->{$row}}) != 0;
		}
		close(OUT);
	}
}

sub write_sparse {
	my $self = shift;
	my $output = shift;

	#use nsort for small key sets ?? Maybe Sort::Key?
	my @split = sort keys %{$self->{countresults}->{seen}->{splitby}};
	my @rows  = sort keys %{$self->{countresults}->{seen}->{rows}};
	my @cols  = sort keys %{$self->{countresults}->{seen}->{cols}};


	#write to file or stdout (if not splitby and no --output)
	foreach my $splitfile (@split) {
		#FIXME  if somebody names it splitfile splitby we're in deep trouble....
		if ($output) {
			my $filename = $output;
			$filename =~ s/(\.txt)?$/-$splitfile.txt/i 
			   unless $splitfile eq "nohit_splitby" && not exists $self->{template}->{$self->{settings}->{split}->{splitby}};
			#add .txt to filename if necessary
			$filename .= "\.txt" unless $filename =~ /\.txt$/i;
			open(OUT,">",$filename);
			select(OUT);
		}
		#print header
		#say join("\t", "name", @cols);
		foreach my $row (@rows) {
			foreach my $key (keys %{$self->{countresults}->{$splitfile}->{$row}}) {
				say join("\t",  $row, $key, $self->{countresults}->{$splitfile}->{$row}->{$key});
			}
		}
		close(OUT);
	}
}

sub print_stats {
	my $self = shift;

	return unless exists $self->{stats};
	print STDERR join("\t", $_, $self->{stats}->{$_}), "\n", foreach keys %{$self->{stats}};
}

1;

