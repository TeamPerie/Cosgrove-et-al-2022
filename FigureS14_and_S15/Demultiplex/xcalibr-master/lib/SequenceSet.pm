use strict;
use warnings;

package SequenceSet;

use Carp;

sub new {
	my $class = shift;
	my $len = shift;
	my $self = {};
	bless $self, $class;

	$self->{seqs} = {};
	$self->{mismatch} = {};
	$self->{len} = $len if defined $len;

	return $self;
}

sub all_seqs {
	my $self = shift;
	return (keys %{$self->seqs}, keys %{$self->mismatch});
}

sub all_names {
	my $self = shift;
	my @names = map { join(":", @$_) } values(%{$self->{seqs}});
	return (@names, map { "mismatch1_".$_ } @names);
}

sub mismatches1 {
	my $self = shift;
	my $sequence = shift;

	my @Nts = split //, $sequence;

	my $ret;
	for my $i (0..$#Nts) {
        for my $mm (map { $_ eq $Nts[$i] ? () : $_ } qw/A T C G N/) {

			my $seq = join("", ($i == 0 ? () : @Nts[0..($i-1)]), $mm, ($i == $#Nts ? () : @Nts[($i+1)..$#{Nts}]));
			#print STDERR "$seq already observed as sequence\n" if exists $self->{seqs}->{$seq};
			#print STDERR "$seq already observed as mismatch\n" if exists $self->{mismatch}->{$seq};
			push @$ret, $seq;
		}
    }
	return $ret;
}

sub read_from_fasta {
	my $self = shift;
	my $file = shift;
	my $mismatches = shift;
	$mismatches = 0 if defined $mismatches;

	croak "File not found $file" unless -e $file;

	require Bio::SeqIO;

	my $in  = Bio::SeqIO->new(-file => $file, -format => 'Fasta');
	while (my $seq = $in->next_seq() ) {
		my $s = $seq->seq;
		$s = substr $s,0,$self->{len} if exists $self->{len};
		unshift @{$self->{seqs}->{$s}}, $seq->id;
		if (defined $mismatches) {
			print STDERR "Sequence $s corresponds to a mismatch\n" if exists $self->{mismatch}->{$s};
			my $r = $self->mismatches1($s);
			$mismatches += @$r;
			push (@{$self->{seqs}->{$_}}, "mismatch1_".$seq->id) foreach @$r;
		}
	}
	print STDERR "Including matching against $mismatches mismatches\n" if defined $mismatches;
}

sub match {
	my $self = shift;
	my $seq = shift;
	if (exists $self->{seqs}->{$seq}) {
		return join(":", @{$self->{seqs}->{$seq}});
	} elsif (exists $self->{mismatch}->{$seq}) {
		return "mismatch1_".join(":", @{$self->{mismatch}->{$seq}});
	}
	return undef;
}

1;

