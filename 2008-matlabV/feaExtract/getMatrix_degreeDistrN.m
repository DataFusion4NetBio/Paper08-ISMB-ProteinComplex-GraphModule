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
% This function gets the nodes' degree distribution to get the network distribution 
% Function Assume the input matrix is the 0/1 sparse graph matrix
% 
% function [ features] = getMatrix_degreeDistrN( maxdegree, matrixgraph, degreePlotTitle , outputPre , opFlag)
% 

function [ features, degree] = getMatrix_degreeDistrN( maxdegree, matrixgraph, degreePlotTitle , outputPre , opFlag)


degree = sum( matrixgraph, 2); 


if ( opFlag )

temp = degree(degree~=0); 
[N, I] = hist(temp, 1:1:maxdegree);
    
% output into figure
graphtitle = sprintf('%s.degree distribute ', degreePlotTitle)
figure
H = loglog(I,N,'.')
xlabel(sprintf('node degree'))
ylabel('cout - frequency ')
title(graphtitle)
errfig = sprintf('%s.degreeDist.fig', outputPre); 
saveas(H, errfig, 'fig'); 



% -------  some statistical numbers ------
disp('In the following, we try to display some statitical summary about the network degree. ')

degree0numEdge = sum(degree == 0); 
disp(' num of nodes having degree 0 '); 
disp(degree0numEdge); 

degree1numEdge = sum(degree == 1); 
disp(' num of nodes having degree 1 '); 
disp(degree1numEdge); 

degree2numEdge = sum(degree == 2); 
disp(' num of nodes having degree 2'); 
disp(degree2numEdge); 

end 


% mean degree 
meandegree = mean(sum(matrixgraph ~=0 )); 
% disp(' mean degree of nodes in average '); 
% disp(meandegree); 

% std of degree 
stddegree = std(full(sum(matrixgraph ~=0 )));
% disp(' std of average degree of nodes'); 
% disp(stddegree); 


% median degree 
mediandegree = median(sum(matrixgraph ~=0 )); 
% disp(' median degree of nodes in average '); 
% disp(mediandegree); 


%max degree 
maxdegree = max(sum(matrixgraph ~=0 )); 



if ( opFlag )
    
% graph density 
graphdensity = nnz(matrixgraph)/prod(size(matrixgraph)); 
disp(' graph density of the sparse matrix '); 
disp(graphdensity); 

% number of edges  
numEdges = nnz(matrixgraph)/2; 
disp(' number of edges of the sparse matrix '); 
disp(numEdges); 

% number of nodes  
numNodes = size(matrixgraph, 1); 
disp(' number of nodes of the sparse matrix '); 
disp(numNodes); 


% ---------------   
% output into file 
outfile = sprintf('%s.degreeDist', outputPre)
fid = fopen(outfile,'w');
fprintf(fid, 'numEdges: %d\n', full(numEdges)); 
fprintf(fid, 'numNodes: %d\n', full(numNodes)); 
fprintf(fid, 'degree0numNodes: %d\n', full(degree0numEdge)); 
fprintf(fid, 'degree1numNodes: %d\n', full(degree1numEdge)); 
fprintf(fid, 'degree2numNodes: %d\n', full(degree2numEdge)); 
fprintf(fid, 'meandegree: %.3f\n', full(meandegree)); 
fprintf(fid, 'stddegree: %.3f\n', full(stddegree)); 
fprintf(fid, 'mediandegree: %d\n', full(mediandegree)); 
fprintf(fid, 'graphdensity: %.6f\n\n', full(graphdensity)); 

partI = find(N ~= 0); 
fprintf(fid,'%d\t%d\n', [ I(partI); N(partI) ]);
fclose(fid);

% --------------
% output related features from this property

mini = 1e-10; 
features = [ full(numNodes); 
             full(numEdges); 
             full(degree1numEdge)/ (full(numNodes) + mini); 
             full(degree2numEdge)/ (full(numNodes) + mini); 
             full(meandegree); 
             full(stddegree); 
             full(mediandegree); 
             full(graphdensity)   ]; 

nodeSize = full(numNodes);
edgeSize = full(numEdges);


else 
    
features = [ full(meandegree); 
             full(stddegree); 
             full(mediandegree); 
             full(maxdegree)   ]; 

end 
