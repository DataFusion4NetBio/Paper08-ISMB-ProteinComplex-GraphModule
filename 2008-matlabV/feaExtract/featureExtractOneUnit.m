%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% wrapper for feature extraction for a small subgraph 
%
% function [ features ] = featureExtractOneUnit( graphMatrix, featureChoiceFlag, graphWeightflag) 
% 
% - featureChoiceFlag {'30fea', '29fea', '29feasvd','density', 'densityNode','33feasvdinfo'}
%

function [ features ] = featureExtractOneUnit( graphMatrix, featureChoiceFlag, graphWeightflag) 


if (strcmp(featureChoiceFlag, 'density') | strcmp(featureChoiceFlag, 'densityNode'))
    [ features ] = generateOneUnitDensity( graphMatrix, graphWeightflag, featureChoiceFlag)  ;     
else % '30fea' or '29fea' or '29feasvd'
    [ features ] = generateOneUnitFeatures( graphMatrix, graphWeightflag, featureChoiceFlag)  ;     
end 