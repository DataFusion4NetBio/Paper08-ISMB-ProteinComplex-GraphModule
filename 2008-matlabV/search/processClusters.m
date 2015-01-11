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
%function [ newClusters ] = processClusters( Clusters , nodeSpaceSize)
%


function [ newClusters ] = processClusters( Clusters , nodeSpaceSize)

global CLUSTER_COLM  CLUSTER_LIMIT  SUBGRPH_SCORE_MIN SUBGRPH_SCORE_MINB  CLSOVERLAPLIMIT; 


numclusters = size(Clusters, 1); 
flag = zeros(numclusters, 1); 

newClusters= cell(0,CLUSTER_COLM);
count = 1; 

if (numclusters == 1)
    newClusters = Clusters; 
    return; 
end 


for i = 1 : (numclusters)
    if ( flag(i) == 0 )    
        curCluster = Clusters{i, 1}; 
        curScore = Clusters{i, 2};

        curV = zeros(nodeSpaceSize, 1); 
        curV(curCluster) = 1; 
        for j = (i + 1) : (numclusters)
          if ( flag(j) == 0 )    
            nexCluster =  Clusters{j,1};
            nexScore = Clusters{j, 2};
           
            nextV  = zeros(nodeSpaceSize, 1); 
            nextV(nexCluster) = 1; 
        
            simil = sum(curV .* nextV)/ min( length(curCluster), length(nexCluster)); 
            overlappednum = sum(curV .* nextV); 
            
            if (( overlappednum >= 2) & ( simil > CLSOVERLAPLIMIT ))
                flag(j) = 1; 
                if  ( flag(i) ~= 1) 
                    flag(i) = 1; 
                    if (length(curCluster) > length(nexCluster))
                        newClusters{count, 1} = curCluster; 
                        newClusters{count, 2} = curScore; 
                        count = count + 1; 
                    else 
                        newClusters{count, 1} = nexCluster; 
                        newClusters{count, 2} = nexScore; 
                        count = count + 1;                     
                   end 
                end 
            end 
          end 
        end 
    end 
    if ( flag(i) == 0 )
        newClusters{count, 1} = curCluster; 
        newClusters{count, 2} = curScore; 
        count = count + 1;             
        flag(i) = 1; 
    end 
end 