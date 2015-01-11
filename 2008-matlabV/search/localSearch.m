%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [  ] = localSearch( graphF , seedChoice, searchParameters, complexScorecutoff,  trainMethodchoice, trainedModelF, featureChoiceFlag, graphWeightflag, clusterOutF , proteinListF )

global CLUSTER_COLM  CLUSTER_LIMIT  SUBGRPH_SCORE_MIN SUBGRPH_SCORE_MINB  CLSOVERLAPLIMIT  CLUSTER_LIMITMIN  CHECKNEIULIM  SIMUT  SIMUALPHA  COMPLEXSMAX; 
CLUSTER_COLM = 2; 
CLUSTER_LIMITMIN = 3; 
SUBGRPH_SCORE_MIN = -100; 
SUBGRPH_SCORE_MINB = -200; 

if (size(searchParameters,1) == 2)
    CLUSTER_LIMIT = 8; 
    CLSOVERLAPLIMIT = 0.7; 
    CHECKNEIULIM = 20; 
    SIMUT = 2
    SIMUALPHA = 2; 
    COMPLEXSMAX = 6; 
elseif  (size(searchParameters,1) ==5 )
    CLUSTER_LIMIT =  searchParameters{3};
    CLSOVERLAPLIMIT =  searchParameters{4};
    CHECKNEIULIM =  searchParameters{5};    
    SIMUT = 2
    SIMUALPHA = 2;     
    COMPLEXSMAX = 6;     
else 
    CLUSTER_LIMIT =  searchParameters{3};
    CLSOVERLAPLIMIT =  searchParameters{4};
    CHECKNEIULIM =  searchParameters{5};      
    SIMUT = searchParameters{6};
    SIMUALPHA = searchParameters{7};    
    COMPLEXSMAX = searchParameters{8};          
end

numSeedtoStart = searchParameters{1};
startnodefile = searchParameters{2};


load(graphF);
graphMatrix = resultingMatrix; 


nodeSize = size(graphMatrix, 1); 
degree = sum(graphMatrix ~= 0);
[sortdegree, sortdegreeI] = sort(degree); 

% -------  Initialized ------------
if ( strcmp(seedChoice, 'randomNode') == 1)
    reOrgNode = randperm(nodeSize); 
    [startingNodes] = reOrgNode(1:numSeedtoStart);
    
elseif ( strcmp(seedChoice, 'allNode') == 1)
    [startingNodes] = sortdegreeI ;    
    
elseif ( strcmp(seedChoice, 'degreeNode') == 1)
    [startingNodes] = sortdegreeI((end-numSeedtoStart+1):end) ;    
    
elseif ( strcmp(seedChoice, 'fileNode') == 1)    

    tempIndexF = '../../temp/temp.curindexf'; 
    logF = '../../temp/temp.log'; 
    FileName = startnodefile; 
    cmd = '!perl ../search/TransferOrf2IndexNew.pl '; 
    command = sprintf('%s  %s  %s  %s > %s ', cmd, proteinListF, FileName, tempIndexF, logF ); 
    eval( command ); 
    oldindex = load( tempIndexF ); 
    Nodes = unique(oldindex(find(oldindex)));    
    
    curNodesDegree = degree(Nodes);
    [sortCurdegree, sortCurdegreeI] = sort(curNodesDegree);
    startingNodes = Nodes(sortCurdegreeI);
end 


disp(sprintf('\n------starting node sets include %d nodes ------------\n', length(startingNodes)));
backupclusteroutF = sprintf('%s.eachCluster', clusterOutF);

% -------  search clusters ------------

count = 0; 

for i = length(startingNodes):-1:1
    curNode = startingNodes(i); 
    disp(sprintf('-- search local cluster for node- %d (positionInFullList: %d)', i, curNode));    
    [ foundFlag ] = findLocalcluster4Node(trainMethodchoice, complexScorecutoff, curNode, graphMatrix, trainedModelF,  featureChoiceFlag, graphWeightflag, backupclusteroutF);
    if (foundFlag)
        count = count + 1; 
        disp(sprintf('-- now we have %d clusters detected.  ', count ))    
    end 
end 


% ---- read in clusters -------
[clusters] = readinEachClusterF(backupclusteroutF);

% ---- remove clusters overlapping -------
updateClusters = processClusters( clusters , nodeSize); 

newClustersNum = size(updateClusters, 1); 
disp(sprintf('\n-- after removing overlapped clusters, now %d clusters ---', newClustersNum))
outputClusters( updateClusters, clusterOutF); 
