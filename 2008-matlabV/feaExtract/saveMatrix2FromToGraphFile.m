function [] = saveMatrix2FromToGraphFile(graphMatrix, SnapGraphfile) 

[i,j,s] = find(graphMatrix); 
proteinSize = size(graphMatrix, 1); 
edgeSize = nnz(graphMatrix); 

% output into file 
fid = fopen(SnapGraphfile, 'w');
fprintf(fid,'# Undirected Graph\n');
fprintf(fid,'# Nodes: %d ; Edges: %d \n', proteinSize, edgeSize );
fprintf(fid,'%d\t%d\n', [ i' ; j' ]);
fclose(fid);
