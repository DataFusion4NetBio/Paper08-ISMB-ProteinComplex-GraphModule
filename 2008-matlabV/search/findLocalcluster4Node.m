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
% function [  flag ]= findLocalcluster4Node(trainMethodchoice, complexScorecutoff, curNode, graphMatrix, trainedModelF,  featureChoiceFlag, graphWeightflag, eachClusterF)
%

function [ flag ]= findLocalcluster4Node(trainMethodchoice, complexScorecutoff, curNode, graphMatrix, trainedModelF,  featureChoiceFlag, graphWeightflag, eachClusterF)

global CLUSTER_COLM  CLUSTER_LIMIT  SUBGRPH_SCORE_MIN CLSOVERLAPLIMIT CLUSTER_LIMITMIN; 

nodeSpaceSize = size(graphMatrix, 1); 

% ---  to avoid to revisit some visted clusters -----


% ---- search local cluster for the node interested -------
curCluster = cell(1, CLUSTER_COLM); 
curCluster{1,1} = [curNode]; 
curCluster{1,2} = SUBGRPH_SCORE_MIN; 

runNum = 0; 
changeflag = 1; 
while  (( runNum < CLUSTER_LIMIT ) & (changeflag))
    [ newCurCluster, changeflag ] = updateOneCluster(runNum, trainMethodchoice, complexScorecutoff, curCluster, graphMatrix, nodeSpaceSize, trainedModelF, featureChoiceFlag, graphWeightflag) ;    
    curCluster =  newCurCluster ; 
    runNum = runNum + 1; 
end 


% -- update the clusters with current found --- 
fid = fopen(eachClusterF, 'a'); 
curNodes = curCluster{1, 1}; 

flag = 0; 
if (length(curNodes) >= CLUSTER_LIMITMIN)
    fprintf(fid,'%d', curNodes(1));    
    for j = 2:length(curNodes)
        fprintf(fid,',%d', curNodes(j));
    end
    fprintf(fid,'  %.5f\n', curCluster{1,2});      
    flag = 1;     
end 
fclose(fid); 
disp(sprintf('-- the most recent cluster size %d  --   ', length(curNodes) ))    
