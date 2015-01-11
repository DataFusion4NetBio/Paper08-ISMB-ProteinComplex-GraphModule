#
# This is a program to generate the average feature related to protein size/length for a input of ORFs 
# 
# Usage: command ProteinsInfoFile OrfFile outFeaFile 
# 


use strict; 
die "Usage: command ProteinsInfoFile OrfFile outFeaFile \n" if scalar(@ARGV) < 3 ;

my ( $ProteinsInfoFile, $OrfListFile, $outFeaFile ) = @ARGV;



#--------------------- read in the ProteinsInfoFile list file -------------------------

my @per_line = (); 
my %proteinsLength = (); 
my %proteinsWeight = (); 
my $count = 0; 

open(PR1, $ProteinsInfoFile) || die(" Can not open file(\"$ProteinsInfoFile\").\n");
while (<PR1>)
{
	chomp; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines

	$count = $count + 1; 
	@per_line = split('\t', $_); 
	$proteinsLength{ "$per_line[0]" } = $per_line[5] ; 
	$proteinsWeight{ "$per_line[0]" } = $per_line[2]/100 ;
}
close(PR1); 
print "# The number of proteins having length infro is : ". keys ( %proteinsLength ) ."\n";
print "# The number of proteins having weight infro is : ". keys ( %proteinsWeight ) ."\n";




#--------------------- read the OrfListFile and process into the features -------------------------

my $weighSum = 0; 
my $lengthSum = 0; 
my $weighMax = 0; 
my $lengthMax = 0; 


open(IN, $OrfListFile) || die(" Can not open file(\"$OrfListFile\").\n");
$count = 0; 
while (<IN>)
{
	chomp; 
	next if /^#/;			#ignore comments
	next if /^$/; 			#ignore blank lines

	$count = $count + 1; 
	my $curorf = $_ ; 
	if ( defined $proteinsLength{ "$curorf" })
	{
		my $cur = $proteinsLength{ "$curorf" }; 
		$lengthSum = $lengthSum + $cur; 
		if ($cur > $lengthMax)
		{	$lengthMax = $cur;  }
	}
	if ( defined $proteinsWeight{ "$curorf" })
	{
		my $cur = $proteinsWeight{ "$curorf" }; 
		$weighSum = $weighSum + $cur; 
		if ($cur > $weighMax)
		{	$weighMax = $cur;  }
	}	
}
close(IN); 
print "Input sugroup: $count ORFs\n";


# ----------------------------------------------
open(OUT, "> $outFeaFile") || die(" Can not open file(\"$outFeaFile\").\n");
#print OUT "#AveLength,MaxLength,AveWeight,MaxWeight\n";
my ($aveL , $aveW);

if ( $count > 0 )
{
	$aveL = $lengthSum/$count; 
	$aveW = $weighSum/$count; 
	print OUT "$aveL,$lengthMax,$aveW,$weighMax\n"; 
}
else {
	print OUT "0,0,0,0\n"; 
}
close(OUT); 
