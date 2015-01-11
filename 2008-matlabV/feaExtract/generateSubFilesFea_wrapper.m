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
% function [] = generateSubFilesFea_wrapper( subgraphsSumF, subgraphEachFileNpre, graphMatrixF,  proteinListF, outFeatureFile, mincomplexSize2Consider, featureChoiceFlag, graphWeightflag  )
%
% - featureChoiceFlag {'30fea', '29fea', '29feasvd','density', 'densityNode','33feasvdinfo'}
%

function [outFeatureFileComplexInf] = generateSubFilesFea_wrapper( subgraphsSumF, subgraphEachFileNpre, graphMatrixF,  proteinListF, outFeatureFile, minSizeconsider, featureChoiceFlag, graphWeightflag  )


% read in the graphMatrix 
data = load(graphMatrixF); 
graphMatrix = data.resultingMatrix; 


% ---- read in the complexes name index --- 
[complex, fileindex, complexsize] = textread(subgraphsSumF, '%s  %d  %d', 'headerlines', 1); 
useComplexIndex = fileindex(find(complexsize >= minSizeconsider)) ; 
NumTrainComplex = length(useComplexIndex) 

trainComplexIndex = useComplexIndex ; 


% ---- Derived the features for the complex units -------  the generated feature file would be comma seperated (CSV) ---
trainPSize = length( trainComplexIndex ); 
fileNames = cell(trainPSize, 1); 
for i =1: trainPSize
    ComplexFileName = sprintf('%s%d', subgraphEachFileNpre, trainComplexIndex(i)); 
    fileNames{i} =  ComplexFileName;    
end 

[numFeature, useFileI ] = generateSubFilesFea( fileNames, graphMatrix,  proteinListF, outFeatureFile, minSizeconsider, featureChoiceFlag, graphWeightflag) ; 


% output the complexes having features into a separate information file 
% complex, fileindex, complexsize   useFileI    outFeatureFile
outFeatureFileComplexInf = sprintf('%s.complexInfo', outFeatureFile); 
outFeatureFileComplexSize = length(useFileI); 

fid = fopen(outFeatureFileComplexInf,'w');
fprintf(fid, '#complexname	complexID	orfsize\n'); 
for i = 1: outFeatureFileComplexSize
    curName = complex{ useFileI(i) };
    curID = fileindex( useFileI(i) ); 
    curSize = complexsize( useFileI(i) ); 
    fprintf(fid,'%s\t%d\t%d\n', curName, curID, curSize);
end 
fclose(fid);