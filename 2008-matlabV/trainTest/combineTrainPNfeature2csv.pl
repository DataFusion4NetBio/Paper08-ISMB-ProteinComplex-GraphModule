# This program is to combine the Feature Pos_file and Neg_file together together into a CSV format feature_lable file
# The input two feature files are assumed to have CSV format
# 

use strict; 
die "Usage: command inPosFeaset inNegFeaset outCSVfile\n" if scalar(@ARGV) < 3;

my ($inPosFeaset, $inNegFeaset, $outCSVfile ) = @ARGV;

my $line_num = 0; 
my $counterp = 0; 	
my $countern = 0; 

open(OUT, "> $outCSVfile") || die(" Can not open file(\"$outCSVfile\").\n"); 

open(IN, $inPosFeaset) || die(" Can not open file(\"$inPosFeaset\").\n"); 
while (<IN>)	
{
	chomp; 
	my $per_line = $_; 
	print OUT $per_line.",1\n"; 
	$counterp = $counterp  + 1; 
	$line_num = $line_num +1; 
}
close(IN); 

open(IN, $inNegFeaset) || die(" Can not open file(\"$inNegFeaset\").\n"); 
while (<IN>)	
{
	chomp; 
	my $per_line = $_; 
	print OUT $per_line.",-1\n"; 
	$countern = $countern  + 1; 
	$line_num = $line_num +1; 
}
close(IN); 
close(OUT); 

print "Total $line_num lines. $counterp pos lines + $countern neg lines.\n"; 