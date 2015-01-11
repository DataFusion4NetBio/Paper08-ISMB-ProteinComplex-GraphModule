%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [clusters] = readinEachClusterF(eachClusterF)


global CLUSTER_COLM ; 
CLUSTER_COLM = 2; 

[nodes, score] = textread(eachClusterF, '%s %f'); 

numclusters = size(score, 1);
clusters = cell(numclusters, CLUSTER_COLM);

for i = 1:numclusters
    clusters{i, 2} = score(i); 
    curArray= zeros(0,1);
    count = 0;
    remain = nodes{i};
    while 1
        [str, newremain] = strtok(remain, ',');
        remain = newremain  ; 
        if isempty(str)
            break;  
        end
        count = count + 1;        
        curArray(count) = str2num(str);
    end
    clusters{i, 1} = curArray; 
end 
