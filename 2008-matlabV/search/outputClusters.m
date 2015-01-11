%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [] = outputClusters( clusters, clusterOutF); 

clusterScoref = sprintf('%s.clusterScore', clusterOutF); 
fidscore = fopen(clusterScoref,'w');
fid = fopen(clusterOutF,'w');

for i =1:size(clusters , 1)
    curCluster = clusters{i, 1}; 
    fprintf(fid,'%d', curCluster(1));
    for j = 2:length(curCluster)
        fprintf(fid,',%d', curCluster(j));
    end
    fprintf(fid,'\n');    
    
    fprintf(fidscore,'cluster%d %d', i, curCluster(1));
    for j = 2:length(curCluster)
        fprintf(fidscore,',%d', curCluster(j));
    end
    tempscore = clusters{i, 2};
    fprintf(fidscore,'  %.5f\n', tempscore);        
end

fclose(fid); 
fclose(fidscore); 