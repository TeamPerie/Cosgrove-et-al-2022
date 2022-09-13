#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: indexhits.pl
#
#        USAGE: ./indexhits.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: YOUR NAME (), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 04/15/2015 01:30:07 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my $index = shift or die "$0 <index>\n";
my $failed = shift;

my $re = $failed ? qr/.*1:[YN]:0:$index$/ : qr/.*1:[N]:0:$index$/;
print STDERR "Including failed reads\n" if $failed;


while(my $line = <>) {
	my $s = <>;
	my $qh = <>;
	my $q = <>;

	if($line =~ /$re/) {
		print $line, $s, $qh, $q;
	}
}

