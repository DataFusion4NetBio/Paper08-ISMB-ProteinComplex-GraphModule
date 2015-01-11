%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



addpath('../feaExtract'); 
addpath('../search'); 
addpath('../trainTest'); 


% --- 0..process your PPI graph and related nodes into two files 

% GraphMatrixF - assumes to be a matlab matrix file including the weighted PPI network  
% MatrixNodeListF - assumes to be a text file with each line including the ORF name 



% ---  1. feature extract for training ----

% minSizeconsider = 3; 
% trainSeperateFileIndexF - assumes to be an index file about all the training complexes and non-complexes
% we assume that each complex or non-complex are in a separate file. The file contains all the related proteins this subgraph relates (ORF name). 
% trainFilePre - where the training subgraph files are 
% featureChoiceFlag - provides multiple choices for generating the subgraph features 

[FeatureInfoFile] = generateSubFilesFea_wrapper( trainSeperateFileIndexF, trainFilePre, GraphMatrixF,  MatrixNodeListF, outFeatureFile, minSizeconsider, featureChoiceFlag  );


% ---  2. training ----


% derivedModelF - output model file name  
% trainMethodchoice {'svm', 'bc'}

complexTrainModel(trainCSVFile, derivedModelF, featureChoiceFlag, trainMethodchoice, trainMethodpara )




% ---  3. searching ----

% seedChoice - which kinds of nodes list to start the searching 

searchParameters = cell(8,1);
searchParameters{1} = 0;
searchParameters{2} = startnodefile ;
searchParameters{3} = 20;   %    CLUSTER_LIMIT 
searchParameters{4} = 0.75; %   CLSOVERLAPLIMIT 
searchParameters{5} = 20; %    CHECKNEIULIM    
searchParameters{6} = 2;   % SIMUT   
searchParameters{7} = 1.8; %    SIMUALPHA    
searchParameters{8} = 2.8; %    COMPLEXSMAX

% clusterOutF - where to ouput the searched clusters

% localSearchParaF  - output the related parameters in this run 
outputParaFile( localSearchParaF,  searchParameters , trainSeperateFileIndexF, trainMethodchoice, trainMethodpara, featureChoiceFlag, seedChoice);
localSearch( GraphMatrixF, seedChoice, searchParameters, traincomplexScorecutoff, trainMethodchoice, derivedModelF, featureChoiceFlag, matrixgraphWeightflag, clusterOutF, MatrixNodeListF); 


