#
# Program to train-test a Bayesian classfier for the 12PPI-graph version 
# graph derived from the DIPS-SVM-detailed
# 
# perl command posClassPriorPara train_CSV_file  model_file nodePosDistFlag[1,0]
#


use strict;
die "Usage: command posClassPriorPara TrainCSVfile feaChoiceFlag model_file nodePosDistFlag[1,0] \n" if scalar(@ARGV) < 5;
my ( $posClassPriorPara, $TrainCSVfile, $feaChoiceFlag, $model_file, $nodePosDistFlag ) = @ARGV;




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




# -------- initialize the parameter array ------------------


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



# -------- read in training file ------------------

my $countp = 0; 
my $countn = 0; 
my $count = 0; 

my ( $curV, $segName, $nodesegName, @data_array); 

open(IN, $TrainCSVfile) || die(" Can not open file(\"$TrainCSVfile\").\n"); 
while (<IN>) 
{
	chomp;
	$count = $count + 1; 
  	@data_array = split(',', $_) ;
	
	if ($#data_array != $numGraphFea )
	{	die "The input 12PPI graph descriptors should have  $numGraphFea + 1 fea version. \n"; 	}
	
	my $flag = $data_array[$#data_array]; 
	if ( $flag == 1 )  
	{	$countp = $countp + 1; 		
		$ci = 1; 
	}
	elsif ( $flag == -1 )  
	{ 	$countn = $countn + 1; 	
		$ci = 0; 	
	}
	else 
	{  	print "Wrong example in CSV train file ! \n"; 	}
	
	# First NodeSize feature 
	$i = 0; 
	$curV = $data_array[$i]; 
	$nodesegName = &getSeqBound($curV, \@nodeNumSeg1 ); 	
	$nodePara{$ci}{"$nodesegName"} ++ ; 
	
	# Then the remaining features 	
	for ($i = 1; $i < $numGraphFea ; $i ++)
	{
		$curV = $data_array[$i]; 
		$segName = &getSeqBound($curV, $featureSegArray[ $i - 1 ] ); 	
		$featureParaArray[ $i - 1 ]->{"$ci.$nodesegName"}{"$segName"} ++ ; 
	}
}
close(IN); 
print "Total POS train : $countp \n"; 
print "Total Rand train : $countn \n";
print "Total train examples : $count \n";  



# -------- output model file  ------------------
# we use the Bayesian beta prior for the parameters

open(OUT, "> $model_file") || die(" Can not open file(\"$model_file\").\n");
my $p = 1e-5; 
my $mvirtual = 1e-1; 
my @prob = (); 

# my $prior = ($countp + 1)/($countn + $countp + 2); 
print OUT "$posClassPriorPara\n";

my ($curHashRef, $curP); 

# We assume the size distribution of the non-complex is the same as the complex units
if ($nodePosDistFlag == 1) 
{
	my $temp = 1; 
	$curHashRef = $nodePara{$temp};
	@prob = &getNBparaCount( $curHashRef , $mvirtual, $p );

	for ($i = 0; $i < 2; $i ++)
	{
		for ($j = 0; $j <= $#nodeNumSeg1 ; $j ++)
		{
			$curP = $curHashRef->{"$nodeNumSeg1[$j]"}; 
			print OUT "$curP," ; 
		}
		$curP = $curHashRef->{"max"}; 
		print OUT "$curP\n" ; 
	}
}
else 
{
# the size distribution of the non-complex is not the same as the complex units	
	for ($i = 0; $i < 2; $i ++)
	{
		$curHashRef = $nodePara{$i}; 		
		@prob = &getNBparaCount( $curHashRef , $mvirtual, $p );

		for ($j = 0; $j <= $#nodeNumSeg1 ; $j ++)
		{
			$curP = $curHashRef->{"$nodeNumSeg1[$j]"}; 
			print OUT "$curP," ; 
		}
		$curP = $curHashRef->{"max"}; 
		print OUT "$curP\n" ; 
	}
}


for ($i = 0; $i < $numGraphFea -1 ; $i ++)
{
	my $curfeaturePara = $featureParaArray[$i]; 
	
	for ($j = 0; $j < $cases ; $j ++)
	{
		$ci = int($j / $nodeCases); 
		$nj = $j - $ci * $nodeCases;
		if ( $nj <= $#nodeNumSeg1 )  
		{	$nsegName = $nodeNumSeg1[$nj] ;
		}
		else 	{	
			$nsegName = "max"; 
		}		
		print OUT ">$ci.$nsegName\n";

		my $curHashRef = $curfeaturePara->{"$ci.$nsegName"} ; 
		@prob = &getNBparaCount( $curHashRef, $mvirtual, $p );

		my $curFeaSegRef = $featureSegArray[$i] ; 
		for ($k = 0; $k <= $#{$curFeaSegRef} ; $k ++)
		{
			$curP = $curHashRef->{"$featureSegArray[$i][$k]"} ; 
			print OUT "$curP," ; 
		}
		$curP = $curHashRef->{"max"}; 
		print OUT "$curP\n" ; 
	}	
}
close(OUT); 



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


# -------- subrountine to get the NB parameter based on the count ------------------

sub getNBparaCount 
{
	my( $countHashRef, $mvirtual, $p ) = @_;
	my $sum = 0; 
	my @probs = (); 
	my ($key, $value); 
	
   	while ( ($key, $value) = each( %$countHashRef ) ) 
   	{
        	$sum = $sum +  $value;
	}	
   	while ( ($key, $value) = each(%$countHashRef) ) 
   	{
        	$countHashRef->{"$key"} = ($value + $mvirtual * $p)/($sum + $mvirtual);
	}	

 	for $key ( sort keys %$countHashRef ) 
 	{
		push(@probs, $countHashRef->{"$key"}); 
    	}
    	return @probs; 
}
