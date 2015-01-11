%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------  get the degree correlation plot ------------------
% The plotted average degree of the immediate neighbors of proteins with a degree k
% 
% function [ cdegree, averageDegreeCorrelate ] = getmatrix_neighborconnectN( maxdegree, degree, matrixgraph, neighboorplotTitle , outputPre , opFlag)
%

function [ cdegree, features ] = getmatrix_neighborconnectN( maxdegree, degree, matrixgraph, neighboorplotTitle , outputPre, opFlag )

proteins = size(matrixgraph, 1);

% in the following arrage, each line contains (k, sum_neighbor_degree, num_neighbor)
%maxdegree = 300; 
cdegree = zeros(maxdegree, 3);
cdegree(:,1) = (1:1:maxdegree)';  
for i = 1: (proteins)
    k = degree(i); 
    if ( k ~= 0)
        partnerI = find(matrixgraph(i, :) ~= 0 ); 
        partnerDegree = degree( partnerI ); 
        cdegree(k, 1) = k; 
        cdegree(k, 2) = cdegree(k, 2) + sum(partnerDegree)/length(partnerDegree) ; 
        cdegree(k, 3) = cdegree(k, 3) + 1;     
    end 
end 

% output the whole average degreeCorrelation value 
averageDegreeCorrelate = sum(cdegree(:, 2))/(sum(cdegree(:, 3)) + 0.00000000000000000000000000000000001) ;


% the k-degree degree correlation vs. degree plot 
for k = 1:maxdegree
    if ( cdegree(k, 3) ~= 0 )
        cdegree(k, 2) =  cdegree(k, 2)/ cdegree(k, 3); 
    end 
end 

stdDegreeCorrelate = std(cdegree(:, 2)) ;
maxDegreeCorrelate = max(cdegree(:, 2)) ;



features = [ averageDegreeCorrelate; 
             stdDegreeCorrelate; 
             maxDegreeCorrelate  ]; 



if ( opFlag )

% output into file 
cdegreeI = find(cdegree(:, 2) ~= 0); 
outfile = sprintf('%s.degreeCorrelate', outputPre)
fid = fopen(outfile,'w');
fprintf(fid, 'averageDegreeCorrelate: %.3f \n', averageDegreeCorrelate); 
fprintf(fid, 'stdKDegreeCorrelate: %.3f \n\n', stdDegreeCorrelate); 
fprintf(fid,'%d\t%.3f\n', [cdegree(cdegreeI, 1)' ;  cdegree(cdegreeI, 2)' ]);
fclose(fid);


% output into a figure
figure 
graphTitle = sprintf('%s - degree correlation plot', neighboorplotTitle); 
H = loglog( cdegree(:, 1), cdegree(:, 2), '.')
xlabel(sprintf(' node degree'))
ylabel(' average degree of the immediate neighbors of the node degree k')
title(graphTitle)
errfig = sprintf('%s.degreeCorrelate.fig', outputPre); 
saveas(H, errfig, 'fig'); 

end 