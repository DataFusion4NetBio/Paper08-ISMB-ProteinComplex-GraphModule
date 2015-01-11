%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [] = outputParaFile( localSearchParaF,  searchParameters , valiChoice, trainMethodchoice, trainMethodpara, featureChoiceFlag, seedChoice)
             

fid = fopen(localSearchParaF ,'a');
fprintf(fid,'\n#---------------------------------------------------\n');
fprintf(fid,'trainvaliChoice\t%s\n', valiChoice);
fprintf(fid,'trainMethodchoice\t%s\n', trainMethodchoice);
fprintf(fid,'trainMethodpara\t%s\n', trainMethodpara);
fprintf(fid,'featureChoiceFlag\t%s\n', featureChoiceFlag);
fprintf(fid,'seedChoice\t%s\n', seedChoice);
fprintf(fid,'StartChoice\t%d\n',searchParameters{1});
fprintf(fid,'startnodefile\t%s\n',searchParameters{2});
fprintf(fid,'CLUSTER_LIMIT\t%d\n',searchParameters{3});
fprintf(fid,'CLSOVERLAPLIMIT\t%.2f\n',searchParameters{4});
fprintf(fid,'CHECKNEIULIM\t%d\n',searchParameters{5});
fprintf(fid,'SIMUTemp\t%.2f\n',searchParameters{6});
fprintf(fid,'SIMUALPHA\t%.2f\n',searchParameters{7});
fprintf(fid,'#---------------------------------------------------\n');
fclose(fid);