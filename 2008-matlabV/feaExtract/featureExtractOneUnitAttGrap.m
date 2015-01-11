%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% wrapper for feature extraction for a small subgraph to generate both graph features and attribute features 
%
% function [ features ] = featureExtractOneUnit( graphMatrix, featureChoiceFlag, graphWeightflag) 
% 
% - featureChoiceFlag {'30fea', '29fea', '29feasvd','density', 'densityNode','33feasvdinfo'}
%

function [ features ] = featureExtractOneUnitAttGrap( graphMatrix, featureChoiceFlag, graphWeightflag, curSubGraphORFlistF) 


            [ graphfeatures ] = featureExtractOneUnit( graphMatrix, featureChoiceFlag, graphWeightflag);      
                        
             if (strcmp(featureChoiceFlag, '33feasvdinfo'))
                tempInfoF = '../../temp/temp.curSub.pinfo'; 
                cmd = '!perl ../feaExtract/YproteinInfo/deriveProInfoFea.pl  ../feaExtract/YproteinInfo/protein_properties.tab  '; 
                command = sprintf('%s  %s  %s  >  %s.log ', cmd, curSubGraphORFlistF, tempInfoF, tempInfoF ) ;
                eval( command ); 
                pinfo = dlmread( tempInfoF, ',' );                   
                features = [ graphfeatures;  pinfo']; 
            else 
                features = graphfeatures; 
            end
