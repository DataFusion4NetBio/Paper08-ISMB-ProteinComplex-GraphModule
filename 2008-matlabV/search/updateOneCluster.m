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
%function [ newCurCluster, changeflag ] = updateOneCluster(round, trainMethodchoice, complexScorecutoff, curCluster, graphMatrix, nodeSpaceSize, trainedModelF, featureChoiceFlag, graphWeightflag) 
% 

function [ newCurCluster, changeflag ] = updateOneCluster(round, trainMethodchoice, complexScorecutoff, curCluster, graphMatrix, nodeSpaceSize, trainedModelF, featureChoiceFlag, graphWeightflag) 

global CLUSTER_COLM  CLUSTER_LIMIT  SUBGRPH_SCORE_MIN SUBGRPH_SCORE_MINB  CLSOVERLAPLIMIT CLUSTER_LIMITMIN  CHECKNEIULIM  SIMUT  SIMUALPHA  COMPLEXSMAX; 

if (strcmp( graphWeightflag, 'svmscore') == 1)
    weightMaxCutoff = 1.0; 
else 
    weightMaxCutoff = 0.2;     
end 


% ---  get the connected nodes for this current cluster
curNodes = curCluster{1,1}; 
poolConnectNode = zeros(nodeSpaceSize,1); 
for i = 1:size(curCluster,1)
    node = curNodes(i); 
    Icon = find( graphMatrix(node, :) );
    poolConnectNode(Icon) = 1; 
end 
poolConnectNode(curNodes) = 0; 
connectedNode = find(poolConnectNode); 
connectedNodeFlagN = zeros(length(connectedNode) , 1); 



% ---  to prefilter the neighbor nodes list -----
simNeigh2Cluster = graphMatrix( connectedNode, curNodes); 
[perNeighMax] = max(simNeigh2Cluster, [], 2); 
tooLooseIndex = find( perNeighMax < weightMaxCutoff); 
connectedNodeFlagN( tooLooseIndex ) = 1; 
[ sortperNeighMax, sortIndex ] = sort(perNeighMax); 
temp = connectedNode; 
connectedNode = temp(sortIndex(end:-1:1)); 




% -- to prevent two nodes clusters -- 
if ((length(curNodes) == 1)&(length(connectedNode) >= 1))
    curNodes = [ curNodes , connectedNode(1) ] ;
    curCluster{1,1} = curNodes; 
end



% ---  to avoid to revisit some visted clusters -----



% ---  evaluate the candidate node to add to the current cluster
predictedOuputF = '../../temp/temp.testoutput'; 
BestN = 0; 
BestS = SUBGRPH_SCORE_MINB;  

checkNeighupLimit = length(connectedNode); 
if (length(connectedNode) > CHECKNEIULIM)
    checkNeighupLimit = CHECKNEIULIM; 
end 
for j = 1: checkNeighupLimit
    curNeigh = connectedNode(j);
    temp = sum(find( curNodes ==  curNeigh)); 
    if (( connectedNodeFlagN(j) == 0 ) & (temp == 0))
            tempCluster = [ curNodes, connectedNode(j) ]; 
            curMatrix = graphMatrix(tempCluster, tempCluster); 
            [ features ] = featureExtractOneUnit( curMatrix, featureChoiceFlag, graphWeightflag);      
            complexTestOneCluster( trainedModelF, features, predictedOuputF, featureChoiceFlag, trainMethodchoice )
            curScore = load(predictedOuputF) ;
            if (curScore > BestS)
                BestN = connectedNode(j); 
                BestS = curScore;
                tempNodes = [curNodes , BestN];
            end 
    end
end 



% --- output the updated cluster --- 
newCurCluster = cell(1, CLUSTER_COLM);
oldScore = curCluster{1,2}; 

if (oldScore >= BestS)  %simulated annealing 
    %prob = exp( BestS - oldScore )/ (SIMUT * power(SIMUALPHA, round))
    prob = exp( (BestS - oldScore ) * ( COMPLEXSMAX/(oldScore + 1e-5)) )/ (SIMUT * power(SIMUALPHA, round))
    if ((rand(1) < prob) & ( BestS > complexScorecutoff))
        disp(sprintf('- In annealing: oldScore- %f ; bestCurS- %f', oldScore, BestS)); 
        newCurCluster{1,1} = [ curNodes , BestN ]; 
        changeflag = 1; 
        newCurCluster{1,2} = BestS;         
    else 
        changeflag = 0; 
        newCurCluster = curCluster; 
    end 
else 
    newCurCluster{1,1} = [ curNodes , BestN ]; 
    changeflag = 1; 
    newCurCluster{1,2} = BestS; 
end 