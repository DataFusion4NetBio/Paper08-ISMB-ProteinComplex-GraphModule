# This program is to convert CSV format to SVMlight format data set 
#!perl -w inputCSV outputSVML

use strict; 
die "Usage: command in_set outputfile\n" if scalar(@ARGV) < 2;
my ($inputfile, $out_file ) = @ARGV;

my $line_num = 0; 
my $counterp = 0; 	
my $countern = 0; 
my $counter = 0; 

open(IN, $inputfile) || die(" Can not open file(\"$inputfile\").\n"); 
open(OUT, "> $out_file") || die(" Can not open file(\"$out_file\").\n"); 

$line_num = 0; 
while (<IN>)	
{
	chomp; 
	my $per_line = $_; 
	my @items = split(',', $per_line); 
	
	my $flag = $items[$#items]; 
	if ($flag != 1) {
		print OUT "-1"; 
		$countern = $countern  + 1; 
	}
	else {
		print OUT "1"; 
		$counterp = $counterp  + 1; 
	}
	
	for($counter = 0 ; $counter < $#items ; $counter++)
	{
		my $temp = $counter + 1; 
		print OUT " $temp:$items[$counter]"; 		
	}
	print OUT "\n"; 
	$line_num = $line_num +1; 
}

print "\n $line_num lines;   pos: $counterp + rand: $countern \n\n"; 
close(IN); 
close(OUT); 


