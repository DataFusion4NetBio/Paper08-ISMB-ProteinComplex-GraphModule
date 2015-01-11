%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------  get the topological coefficient plot ------------------
% 
% Defined as TCp= average(J(p,j)/kp), with J(p, j) denoting the number of nodes to which both p and j are linked,
% plus 1 if there is a direct link between p and j. J(p,j) is defined for all proteins p (with kp>1) and j in the network
% which share at least one common interaction partner. kp is the number of links of node p. 
% 
% The topological coefficient is thus calculated only for proteins with more than one link  
% 
% TCp is a relative measure for the extent to which a protein shares interaction partners with other proteins 
% in the network and also reflects the number of rectangles that pass through a node.
% 
% The maximum topological coefficient [TCpm= max(J(p,j)/kp)] for every given protein p indicates the protein j with
% the maximum number of shared interaction partners in the network. Pairs of proteins, which share a maximum
% number of interaction partners are thus identified by means of TCpm. A high TCpm value points towards a functional
% connection between the two proteins, although they may not be directly linked.
% 
% The above graph is also in log-log scale 
%  
% function [ topologicc,  features ] = getMatrix_topolcoefficient( maxdegree, degree, matrixgraph, graphtitle, outputPre , opFlag) 
% 

function [ topologicc,  features ] = getMatrix_topolcoefficient( maxdegree, degree, matrixgraph, graphtitle, outputPre , opFlag) 


proteins = size(matrixgraph, 1) ;

cjpj = zeros(proteins, 2);         
for i = 1: proteins
    k = degree(i); 
    if ( k > 1 )
        % find the 1st level neighbor        
        partnerI = find(matrixgraph(i, :) ~= 0 ); 
        
        % find the 2nd level neighbor
        jpj = zeros(proteins, 2);         
        for j = 1: proteins
            if ( j ~= i )
                otherpartnerI =  find(matrixgraph(j, partnerI) == 1); 
                matchedI = partnerI( otherpartnerI ); 
                if ( length(matchedI) > 0)
                    jpj(j, 1) = length(matchedI) + matrixgraph(i, j); 
                    jpj(j, 2) = 1; 
                end 
            end 
        end       
        cjpj(i, 1) = k ;            
        cjpj(i, 2) = sum(jpj(:, 1))/(sum(jpj(:, 2)) * k + 0.000000000000000000000001 ); 
        if (cjpj(i, 2) > 1)
            disp('error - topologicalC large than 1. ')
            keyboard
        end 
    end    
end

topologicc = zeros(maxdegree, 2);
for k = 1 : maxdegree
    ccI = find( cjpj(:, 1) == k ); 
    if (length(ccI) > 0)
        kcc = cjpj(ccI, 2); 
        topologicc(k, 1) = k; 
        topologicc(k, 2) = mean(kcc); 
    end 
end 

averagetopologicc = sum(cjpj(:, 2))/proteins ;
stdtopologicc = std(cjpj(:, 2));
maxtopologicc = max(cjpj(:, 2));

features = [ averagetopologicc ;  stdtopologicc; maxtopologicc ]; 




if (opFlag)

% output into file 
tempI = find(topologicc(:, 2) ~= 0); 
outfile = sprintf('%s.topologicCoefficient', outputPre)
fid = fopen(outfile,'w');
fprintf(fid, 'averagetopologicc: %.3f \n', averagetopologicc ); 
fprintf(fid, 'stdKtopologicc: %.3f \n\n', stdtopologicc ); 
fprintf(fid,'%d\t%.3f\n', [ topologicc(tempI, 1)' ; topologicc(tempI, 2)' ]);
fclose(fid);


% output to figures 
figure 
tctitle = sprintf('%s.topologicalCoefficient ',graphtitle)

H = loglog( topologicc(:, 1), topologicc(:, 2), '.')
xlabel(sprintf(' node degree '))
ylabel(' the average topological coefficient of all proteins with degree k ')
title(sprintf('%s. average', tctitle))
errfig = sprintf('%s.aveTopologicCoefficient.fig', outputPre); 
saveas(H, errfig, 'fig'); 


figure 
H = loglog( cjpj(:, 1), cjpj(:, 2), '.')
hold on 
plot( topologicc(:, 1), topologicc(:, 2), 'ro')
xlabel(sprintf(' node degree '))
ylabel(' the full topological coefficient of all proteins with degree k ')
title(sprintf('%s. whole', tctitle))
errfig = sprintf('%s.fullTopologicCoefficient.fig', outputPre); 
saveas(H, errfig, 'fig'); 

end 