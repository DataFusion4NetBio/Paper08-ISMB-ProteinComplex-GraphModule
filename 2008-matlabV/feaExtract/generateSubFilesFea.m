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
% function [numFeature] = generateSubFilesFea( subgraphFileNames, graphMatrix,  proteinListF, outFeatureFile, featureChoiceFlag, graphWeightflag) 
% 
% - featureChoiceFlag {'30fea', '29fea', '29feasvd','density', 'densityNode','33feasvdinfo'}


function [numFeature, useFileI ] = generateSubFilesFea( subgraphFileNames, graphMatrix,  proteinListF, outFeatureFile, minSizeconsider, featureChoiceFlag, graphWeightflag) 

tempIndexF = '../../temp/temp.curSuggraph.index'; 
templogF = '../../temp/temp.log'; 

fileNames = subgraphFileNames; 
fileNum = size( subgraphFileNames, 1); 
useFileI = zeros(0,1); 
count = 0; 

fid = fopen(outFeatureFile,'w');
for i = 1: fileNum
    cmd = '!perl ../feaExtract/TransferOrf2Index.pl '; 
    command = sprintf('%s  %s  %s  %s >  %s ', cmd, proteinListF, fileNames{i}, tempIndexF, templogF ) ;
    eval( command ); 
    oldindex = load( tempIndexF ); 
  
    index = oldindex(oldindex > 0);     
        
    if (length(index) >= minSizeconsider ) 
        count = count + 1; 
        useFileI(count) = i; 
        
        curMatrix = graphMatrix(index, index); 
    
        [ vfeature ] = featureExtractOneUnitAttGrap( curMatrix, featureChoiceFlag, graphWeightflag, fileNames{i});
        fprintf(fid,'%.6f', vfeature(1) );
        for j =2: length( vfeature )
            fprintf(fid,',%.6f', vfeature(j) );    
        end                 
        fprintf(fid,'\n'); 
        
    end 
end 
fclose(fid);

numFeature = length(vfeature) 