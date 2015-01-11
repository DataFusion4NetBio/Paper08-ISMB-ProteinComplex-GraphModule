%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eigenVfeatures] = getMatrix_svds( graphMatrix ) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate svd 

% number_of_edges = nnz(A);    %This is the total number of edges

%% We could use MATLAB's svd routines. This is slow 

K = 3; 
[ eigenVfeatures ] = svds(graphMatrix, K );


%% could also use the Lanczos algorithm, which is *much* faster.
%% For example, on Epinions:
%%           Number_Eigenvalues   Time(Lanczos)    Time(MATLAB's SVD)
%%                  10                5.5                 26.01
%%                  40               34.7                169.55
%%
%%
% %% need to addpath to the package path before hand: addpath propack;
% [U,S,V] = lansvd(graphMatrix, firstKSingularValue, 'L');   % 'L' says find the largest eigenvalues
% 
% v1o = abs(V(:,1));
% u1o = abs(U(:,1));
% 
% % %if sum(v1o(1:100)) < 0      %Note: In practice, this works all the time
% % %   v1o = -v1o;              %      to make sure elements of v1 and u1
% % %   u1o = -u1o;              %      are positive (except noise with very
% % %end                         %      small values)
% % 
% % lambdas = diag(S); 
% % eigenVfeatures = lambdas; 
% % 
