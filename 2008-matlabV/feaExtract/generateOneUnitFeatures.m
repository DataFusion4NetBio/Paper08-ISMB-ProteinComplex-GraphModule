%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% function [ vfeature ] = generateOneUnitFeatures( graphMatrix, graphWeightflag) 
%
% - featureChoiceFlag {'30fea', '29fea', '29feasvd','density', 'densityNode','33feasvdinfo'}


function [ vfeature ] = generateOneUnitFeatures( graphMatrix, graphWeightflag, featurechoice) 

if ( nargin < 3 )
    featurechoice = '29feasvd';
end 


% 1. 
nodeSize = size(graphMatrix, 1); 
 
if ( nodeSize <=1 ) 
    disp('Bug: subgraph has less than two nodes. '); 
    keyboard;
end 


% 1. 
%graphdensity = nnz(graphMatrix)/prod(size(graphMatrix)); 
graphdensity = nnz(graphMatrix)/(nodeSize*nodeSize - nodeSize); 

% 2. 
meanEdge = sum(sum(graphMatrix))./(nodeSize*nodeSize - nodeSize); 

% 3.
stdEdge = std(full(graphMatrix(:))); 

% here we assume that the edge weight must be positive 
nnSelect = graphMatrix(find( graphMatrix > 0 )); 
if (length(nnSelect) > 0)
    % 4.    
    EdgeMean = mean( full( nnSelect )); 

    % 5.
    EdgeStd = std( full(nnSelect )); 
else 
    EdgeMean = 0; 
    EdgeStd = 0;     
end 

% 6
cutoff = [ 1.0, 1.2,  1.5, 1.8, 2.2, 2.6, 3.0 ] ; 
useCut = cutoff; 

fraction = zeros( length( useCut ) , 1); 
for i = 1:length( useCut )
     nn =  find( graphMatrix > useCut(i)); 
     fraction(i) = size(nn,1)/(nodeSize * nodeSize - nodeSize); 
end 


%7 features from other graph properties 
binGraph = graphMatrix; 
binGraph(binGraph ~= 0) = 1; 

% - featureChoiceFlag {'30fea', '29fea', '29feasvd','density'}
if (strcmp(featurechoice, '29fea') ==1)
    choosePropertyFlag = [1,1,1,1,0,1,0];     
elseif (strcmp(featurechoice, '29feasvd') ==1)
    choosePropertyFlag = [1,1,1,1,0,0,1]; 
elseif (strcmp(featurechoice, '30fea') ==1)
    choosePropertyFlag = [1,1,1,1,1,0,0];     
end 

maxdegree = 30; 
[ features ] = getMatrix_graphPropertyFeaN( binGraph, maxdegree, choosePropertyFlag, 0, '', '' ); 


% construct the feature vector 
wholefeature = [ nodeSize; 
                 graphdensity; 
                 meanEdge; 
                stdEdge; 
                EdgeMean; 
                EdgeStd; 
                fraction;
                features  ]; 

vfeature = full(wholefeature); 

