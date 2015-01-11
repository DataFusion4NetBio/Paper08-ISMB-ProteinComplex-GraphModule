# 
# Program to test a Bayesian classfier for the 12PPI-graph-feature version 
# graph derived from the DIPS-SVM-detailed
# 

use strict;
die "Usage: command  model_file  feaChoiceFlag TestCSVfile   outputfile \n" if scalar(@ARGV) < 4;
my (  $model_file,  $feaChoiceFlag, $TestCSVfile,  $outputfile ) = @ARGV;


# -------- initialize the discretization boundary of each feature ------------------
my @nodeNumSeg1 =( 3, 4, 5, 7, 10 ) ;
my @densitySeg2 =( 0.1, 0.3, 0.5, 0.7, 0.9 ) ;
my @edgeSeg3_6 =( 1, 3, 5 );
my @fraSeg7_13 = ( 0.1, 0.3, 0.5, 0.7, 0.9 ) ;
my @degreSeg14_17 =( 2, 6, 10 ) ;
my @degreCSeg18_20 =( 2, 5, 10 ) ;
my @clusterCSeg21_23 =( 0.39, 0.88 ) ;
my @topCSeg24_26 =( 0.39, 0.88 ) ;
my @eigSeg27 =( 2, 5, 8 ) ;
my @eigSeg28_29 =( 1, 2.5 ) ;
#my @diamSeg27 =( 2, 3, 4 ) ;


my ($i, $j, $k); 
my @featureSegArray = (); 
my $numGraphFea; 

if (($feaChoiceFlag eq '29fea' )||( $feaChoiceFlag eq '29feasvd' )) 
{
	$numGraphFea = 29; 	
	$featureSegArray[0] = [ @densitySeg2 ]; 
	for ($i = 1; $i <= 4; $i++ ) {
		$featureSegArray[ $i ] = [ @edgeSeg3_6 ]; 
	}
	for ($i = 5; $i <= 11; $i++ ) {
		$featureSegArray[ $i ] = [ @fraSeg7_13 ]; 
	}
	for ($i = 12; $i <= 15; $i++ ) {
		$featureSegArray[ $i ] = [ @degreSeg14_17 ]; 
	}
	for ($i = 16; $i <= 18; $i++ ) {
		$featureSegArray[ $i ] = [ @degreCSeg18_20 ]; 
	}
	for ($i = 19; $i <= 22; $i++ ) {
		$featureSegArray[ $i ] = [ @clusterCSeg21_23 ]; 
	}
	for ($i = 22; $i <= 24; $i++ ) {
		$featureSegArray[ $i ] = [ @topCSeg24_26 ]; 
	}
	for ($i = 25; $i <= 25; $i++ ) {
		$featureSegArray[ $i ] = [ @eigSeg27 ]; 
	}
	for ($i = 26; $i <= 27; $i++ ) {
		$featureSegArray[ $i ] = [ @eigSeg28_29 ]; 
	}	
}
elsif ($feaChoiceFlag eq 'densityNode') {
	$numGraphFea = 2; 	
	$featureSegArray[0] = [ @densitySeg2 ]; 		
}
else {
	print "Wrong input parameter for BN classifier: feaChoiceFlag - $feaChoiceFlag \n";
}




# -------  read in the model file  --------------


my %nodePara = (); 
for ($i = 0; $i < 2; $i ++)
{
	my %curPara = (); 
	for ($j = 0; $j <= $#nodeNumSeg1 ; $j ++)
	{
		$curPara{"$nodeNumSeg1[$j]"} = 0; 
	}
	$curPara{"max"} = 0;
	$nodePara{$i} = \%curPara; 
}


my ($ci, $nj, $nsegName ); 
my @featureParaArray = (); 

my $nodeCases = $#nodeNumSeg1 + 2; 
my $cases = 2 * $nodeCases ; 

for ($i = 0; $i < $numGraphFea -1 ; $i ++)
{
	my $curFeaSegRef = $featureSegArray[$i] ; 

	my %featurePara = (); 
	for ($j = 0; $j < $cases ; $j ++)
	{
		my %curPara = (); 
		for ($k = 0; $k <= $#{$curFeaSegRef} ; $k ++)
		{
			$curPara{"$featureSegArray[$i][$k]"} = 0; 
		}
		$curPara{"max"} = 0;
		
		# here $j represent the index outof the items from (C_i, N_j)
		$ci = int($j / $nodeCases); 
		$nj = $j - $ci * $nodeCases;
		if ( $nj <= $#nodeNumSeg1 )  
		{	$nsegName = $nodeNumSeg1[$nj] ;
		}
		else 	{	
			$nsegName = "max"; 
		}		
		$featurePara{"$ci.$nsegName"} = \%curPara; 
	}	
	$featureParaArray[$i] = \%featurePara; 
}








# ---------------------------  begin reading in the parameters
my ( $curline, @items, $count  ) ; 
open(MD, $model_file) || die(" Can not open file(\"$model_file\").\n"); 

$count = 0; 
$curline = <MD>;
chomp($curline); 
$count = $count + 1; 
my $prior = $curline; 
my $curHashRef; 

for ($i = 0; $i < 2; $i ++)
{
	@items = (); 
	$curline = <MD>;
	$count = $count + 1; 
	chomp($curline); 
	@items = split(/,/, $curline); 
	
	$curHashRef = $nodePara{$i}; 	
	for ($j = 0; $j <= $#nodeNumSeg1  ; $j ++)
	{
		$curHashRef->{"$nodeNumSeg1[$j]"} = $items[$j]; 
	}
	$curHashRef->{"max"} =  $items[$#nodeNumSeg1  + 1]; 
}

my $featurePara; 
for ($i = 0; $i < $numGraphFea - 1; $i ++)
{
	$featurePara = $featureParaArray[$i]; 
	
	for ($j = 0; $j < $cases ; $j ++)
	{
		$ci = int($j / $nodeCases); 
		$nj = $j - $ci * $nodeCases;
		if ( $nj <= $#nodeNumSeg1  )  
		{	$nsegName = $nodeNumSeg1[$nj] ;
		}
		else 	{	
			$nsegName = "max"; 
		}		
		
		$curline = <MD>;
		chomp($curline); 
		$count = $count + 1; 
		if ($curline ne ">$ci.$nsegName")
		{	print "A mistake happens at line: $count \n"; } 
	
		$curline = <MD>;
		chomp($curline); 
		$count = $count + 1; 
		@items = split(/,/, $curline); 	

		$curHashRef = $featurePara->{"$ci.$nsegName"} ; 
		my $curFeaSegRef = $featureSegArray[$i] ; 
		for ($k = 0; $k <= $#{$curFeaSegRef} ; $k ++)
		{
			$curHashRef->{"$featureSegArray[$i][$k]"} = $items[$k]; 
		}
		$curHashRef->{"max"} = $items[ $#{$curFeaSegRef} + 1 ]; 
	}	
}
close(MD); 
print "Total $count lines in modle files ! \n"; 




# -------- read in the test file  and predict ------------------

open(IN, $TestCSVfile) || die(" Can not open file(\"$TestCSVfile\").\n"); 
open(OUT, "> $outputfile") || die(" Can not open file(\"$outputfile\").\n");

my ( $curV, $segName, $nodesegName, @data_array); 

my $count = 0; 
my $curP = 0; 
while(<IN>)
{
	chomp;
  	@data_array = split(',', $_) ;
	if ($#data_array != $numGraphFea  )
	{	die "The input 12PPI graph descriptors should have  $numGraphFea fea version. \n"; 	}
	my $flag = $data_array[$#data_array]; 
	$count = $count + 1; 

	my $sumPos = 0; 
	my $sumRand = 0; 
	
	# First NodeSize feature 
	$i = 0; 
	$curV = $data_array[$i]; 
	$nodesegName = &getSeqBound($curV, \@nodeNumSeg1  ); 	

	$ci = 1; 
	$curP = $nodePara{$ci}{"$nodesegName"}; 
	$sumPos = $sumPos + log($curP); 
	$ci = 0; 
	$curP = $nodePara{$ci}{"$nodesegName"}; 
	$sumRand = $sumRand + log($curP); 

	# Then the remaining features 	
	for ($i = 1; $i < $numGraphFea ; $i ++)
	{
		$curV = $data_array[$i]; 
		$segName = &getSeqBound($curV, $featureSegArray[ $i - 1 ] ); 	
		
		$ci = 1; 
		$curP = $featureParaArray[$i -1]->{"$ci.$nodesegName"}{"$segName"}; 
		$sumPos = $sumPos + log($curP); 
		$ci = 0; 
		$curP = $featureParaArray[$i -1]->{"$ci.$nodesegName"}{"$segName"}; 
		$sumRand = $sumRand + log($curP); 		
	}

	$sumPos = $sumPos + log($prior); 
	$sumRand = $sumRand + log( 1- $prior); 	
	#my $predict = exp($sumPos )/( exp($sumRand ) ); 
	my $predict = ( $sumPos - $sumRand ) / 100; 
	print OUT "$predict\n"; 
}

close(IN); 
close(OUT); 

print "Total Test : $count \n"; 



# -------- subrountine to get the SeqBoundarayName based on a value and the range_seqmentation for that feature ------------------

sub getSeqBound 
{
	my($value, $segArraryRef) = @_;
	my @seArray = @$segArraryRef; 
	
	my ($array_element, $i); 
	if ( $value > $seArray[$#seArray] )
	{	return 'max' }; 
	for ($i = 0; $i <= $#seArray; $i ++)
	{
		if ( $value <= $seArray[$i] )
		{	
			return "$seArray[$i]"; 
			last; 
		}
	}
}

