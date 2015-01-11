%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------  get the clustering coefficient plot ------------------
% calculated for every protein that has more than one interaction partner: CCp = 2n/kp(kp-1), 
% with n as the number of links connecting the kp neighbors of node p to each other. 
% 
% C(k) is defined as the average clustering coefficient of all proteins with k links and 
% the negative gradient on a log-log plot indicates hierarchical modularity. 
% 
% function [ clusterc , features ] = getmatrix_clusteringcoefficientN( maxdegree, degree, matrixgraph , ccplotTitle , outputPre, opFlag )
% 

function [ clusterc , features ] = getMatrix_clusteringcoefN( maxdegree, degree, matrixgraph , ccplotTitle , outputPre , opFlag)

proteins = size(matrixgraph, 1);

ckp = zeros(proteins, 2); 
for i = 1: proteins
    k = degree(i); 
    if ( k > 1 )
        partnerI = find(matrixgraph(i, :) ~= 0 ); 
        sumatrix = matrixgraph(partnerI, partnerI);      
        kp = length(partnerI); 
        kn = sum(sumatrix(:)); 
        ccp =  kn /(kp * (kp -1)); 
        ckp(i, 1) = k;
        ckp(i, 2) = ccp;        
    end                    
end

averageClusterc = full(sum(ckp(:, 2))/proteins) ;
stdClusterc = std( ckp(:, 2) ) ;
maxClusterc = max( ckp(:, 2) ) ;

features = [ averageClusterc;  stdClusterc ; maxClusterc]; 




clusterc = zeros(maxdegree, 2);
for k = 1 : maxdegree
    ccI = find( ckp(:, 1) == k ); 
    if (length(ccI) > 0)
        kcc = ckp(ccI, 2); 
        clusterc(k, 1) = k; 
        clusterc(k, 2) = mean(kcc); 
    end 
end 


if (opFlag)


% output into file 
cluscI = find(clusterc(:, 2) ~= 0); 

outfile = sprintf('%s.clusterCoefficient', outputPre)
fid = fopen(outfile,'w');
fprintf(fid, 'averageClusterc: %.3f \n', averageClusterc); 
fprintf(fid, 'stdKClusterc: %.3f \n\n', stdClusterc); 
fprintf(fid,'%d\t%.3f\n', [clusterc(cluscI, 1)' ; clusterc(cluscI, 2)' ]);
fclose(fid);




% output to the plots 
graphtitle = sprintf('%s.clusteringCoefficient ', ccplotTitle)

figure 
H = loglog( clusterc(:, 1), clusterc(:, 2), '.')
xlabel(sprintf(' node degree '))
ylabel(' the average clustering coefficient of all proteins with degree k ')
title(sprintf('%s. average ', graphtitle))
errfig = sprintf('%s.AveClusterCoefficient.fig', outputPre); 
saveas(H, errfig, 'fig'); 


figure 
H = loglog( ckp(:, 1), ckp(:, 2), '.')
hold on 
plot( clusterc(:, 1), clusterc(:, 2), 'ro')
xlabel(sprintf(' node degree '))
ylabel(' the full clustering coefficient of all proteins with degree k ')
title(sprintf('%s. whole', graphtitle))
errfig = sprintf('%s.FullClusterCoefficient.fig', outputPre); 
saveas(H, errfig, 'fig'); 

end 