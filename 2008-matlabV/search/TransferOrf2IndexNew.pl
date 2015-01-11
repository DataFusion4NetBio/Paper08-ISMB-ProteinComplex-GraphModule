#
# This is a program to transform the ORF name into the indexNum based on a predefied ORF list 
# 
# Usage: command ProteinsListFile OrfFile indexFile 
# 


use strict; 
die "Usage: command ProteinsListFile OrfFile indexFile \n" if scalar(@ARGV) < 3 ;

my ( $ProteinsListFile, $OrfFile, $indexFile ) = @ARGV;


#--------------------- read in the proteinName list file -------------------------

my @per_line = (); 
my %proteinsName = (); 
my $count = 0; 

open(PR1, $ProteinsListFile) || die(" Can not open file(\"$ProteinsListFile\").\n");
while (<PR1>)
{
	chomp; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines

	$count = $count + 1; 
	@per_line = split('\t', $_); 
	$proteinsName{ "$per_line[0]" } = $count ; 
}
close(PR1); 
print "# The size of proteins space this partition is on : ". keys ( %proteinsName ) ."\n";


#--------------------- process the ORFfile into the indexfile -------------------------

open(OUT, "> $indexFile") || die(" Can not open file(\"$indexFile\").\n");
open(IN, $OrfFile) || die(" Can not open file(\"$OrfFile\").\n");
$count = 0; 
while (<IN>)
{
	chomp; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines

	$count = $count + 1; 
	@per_line = split('\t', $_); 
	my $curorf = $per_line[0] ; 
	if ( defined $proteinsName{ "$curorf" })
	{
		my $curIndex = $proteinsName{ "$curorf" }; 
		print OUT "$curIndex\n"; 
	}
	else {
		print OUT "0\n"; 		
	}	
}
close(IN); 
close(OUT); 
