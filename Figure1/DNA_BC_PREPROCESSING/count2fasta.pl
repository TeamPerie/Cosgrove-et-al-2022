#!/usr/bin/perl 
#===============================================================================
#
#         FILE: count2fasta.pl
#
#        USAGE: ./count2fasta.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Arno Velds (), a.velds@nki.nl
#      COMPANY: NKI
#      VERSION: 1.0
#      CREATED: 11/09/2012 11:17:00 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

#skip header
my $line = <>;
while(<>){
my @e = split /\t/;
print join("\n", ">$e[0]", $e[0]), "\n";
}
