#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  cuttags.pl
#
#        USAGE:  ./cuttags.pl  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Arno Velds (AV), a.velds@nki.nl
#      COMPANY:  The Netherlands Cancer Institute
#      VERSION:  1.0
#      CREATED:  02/22/2012 08:42:07 AM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;


my $countsfile = shift;
my $blastfile = shift;

open(B,"<", $blastfile) or die "Blastfile!";
open(C,"<",$countsfile) or die "Counts file\n";

my %blast;
while(my $line = <B>) {
	my @e = split /\t/, $line;
	#check if the blast hit is in the forward direction
	next unless $e[9] > $e[8];
	if(exists $blast{$e[0]}) {
		warn "Duplicate blast hit for $e[0], skipping\n";
		print STDERR $line;
		next;
	}
	$blast{$e[0]} = \@e;
}
close(B);

#parse counts file and clip sequences
while(my $line = <C>) {
	my @e = split /\t/, $line;
	if(exists $blast{join("Y",$e[0],$e[1])}) {
		my $new = substr $e[0], 0, $blast{join("Y",$e[0],$e[1])}->[6] -1; 
		print join("\t", $new, $line);
	} else {
		print join("\t", $e[0], $line);
	}
}
close(C);

