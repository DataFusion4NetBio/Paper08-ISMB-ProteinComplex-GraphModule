%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ vfeature ] = generateOneUnitDensity( graphMatrix, graphWeightflag, featureChoiceFlag) 

nodeSize = size(graphMatrix, 1); 
if ( nodeSize <=1 ) 
    disp('Bug: subgraph has less than two nodes. '); 
    keyboard;
end 
    
graphdensity = nnz(graphMatrix)/(nodeSize*nodeSize - nodeSize); 

% construct the feature vector 
if (strcmp(featureChoiceFlag, 'densityNode'))
    wholefeature = [ nodeSize; graphdensity ]; 
else 
    wholefeature = [ graphdensity ];     
end 
vfeature = full(wholefeature); 