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
% function [ nodeSize, edgeSize, features ] = getMatrix_graphPropertyFeaN( graphMatrix, maxdegree, choosePropertyFlag, opFlag, graphtitle , outputPre )
% 
% This function is a wrapper to various graph_property functions (version 2)
% 
% Note: the graph properties we process here works only on the binary graph matrix 
% 
% we choose the graph properties the user request to calculate - based on parameter: 
% Related INPUT: "choosePropertyFlag"
% 
% PropertyFlag: 
% 1- degree distribution 
% 2- degree correlation
% 3- clustering coefficient
% 4- topological coefficient
% 5- eigenAnalysis
% 6- edge stress 
% 7- Snap generatd properties 


function [ features ] = getMatrix_graphPropertyFeaN( graphMatrix, maxdegree, choosePropertyFlag, opFlag, graphtitle , outputPre )


if (choosePropertyFlag(1) == 1 )
    [  degreeFeatures, degree ] = getMatrix_degreeDistrN( maxdegree, graphMatrix, graphtitle, outputPre, opFlag); 
    %keyboard;
    
    if (choosePropertyFlag(2) == 1 )  
        [ cdegree, DegreeCorrelateFeatures ] = getMatrix_neighborconnectN( maxdegree, degree, graphMatrix,  graphtitle , outputPre, opFlag  ); 
        %keyboard;
    end 

    if (choosePropertyFlag(3) == 1 )
        [ clusterc , ClustercFeatures ] = getMatrix_clusteringcoefN(  maxdegree, degree, graphMatrix, graphtitle , outputPre, opFlag);
        %keyboard;
    end 

    if (choosePropertyFlag(4) == 1 )
        [ topologicc,  topologiccFeatures ] = getMatrix_topolcoefficientN( maxdegree, degree, graphMatrix, graphtitle, outputPre, opFlag); 
        %keyboard;
    end     
end 


effectivediameterAvg = zeros(0, 1);
eigenVfeatures = zeros(0, 1);
if (choosePropertyFlag(5) == 1 )
    [ effectivediameterAvg, eigenVfeatures ] = getMatrix_snapPropertyN( graphMatrix, graphtitle , outputPre, opFlag ) ; 
    %keyboard
end 


if (choosePropertyFlag(6) == 1 )
    [ eigenVfeatures ] = getMatrix_snapEig( graphMatrix, graphtitle , outputPre, opFlag ) ; 
    %keyboard
end 

if (choosePropertyFlag(7) == 1 )
    [ eigenVfeatures ] = getMatrix_svds( graphMatrix ) ; 
    %keyboard
end 

features = [    degreeFeatures; 
                DegreeCorrelateFeatures; 
                ClustercFeatures; 
                topologiccFeatures; 
                effectivediameterAvg;
                eigenVfeatures]; 