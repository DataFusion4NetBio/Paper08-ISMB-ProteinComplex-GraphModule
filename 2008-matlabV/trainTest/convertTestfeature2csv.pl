#
# This program is to convert the test Feature into a CSV format feature_lable file
# The input feature files is assumed to have CSV format
# 

use strict; 
die "Usage: command inTestFeaset outCSVfile\n" if scalar(@ARGV) < 2;

my ($inTestFeaset, $outCSVfile ) = @ARGV;

my $line_num = 0; 

open(OUT, "> $outCSVfile") || die(" Can not open file(\"$outCSVfile\").\n"); 
open(IN, $inTestFeaset) || die(" Can not open file(\"$inTestFeaset\").\n"); 
while (<IN>)	
{
	chomp; 
	my $per_line = $_; 
	print OUT $per_line.",1\n"; 
	$line_num = $line_num +1; 
}
close(IN); 
close(OUT); 

print "Total $line_num lines. \n"; 