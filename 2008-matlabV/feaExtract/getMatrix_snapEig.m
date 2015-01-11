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
% function [ eigenVfeatures ] = getMatrix_snapEig( graphMatrix, graphtitle , outputPre , opFlag)  
% 
% This is a matlab wrapper to use the Snap to get the eigenvalue graph properties
% Suppose that the ./Snap directory exists in the local directory 


function [ eigenVfeatures ] = getMatrix_snapEig( graphMatrix, graphtitle , outputPre, opFlag)  

% We first need to transform the graphMatrix into the Snap input file format 
% The format is the (from-to) format: Graph is in edge list format 
% 
% Typical usage of Snap:
%   SnapApp.exe -i:graph.txt -o:OUT -t:"PLOT TITLE" -p:cdhw
% 
% If edges are  undirected (directed) use command line parameter
%   -d:F (-d:T)
% To plot only certain plots use parameter -p:<string> where string is composed of the following letters:
%   c: cummulative degree distribution
%   d: degree distribution
%   h: hop plot
%   w: distribution of weakly connected components
%   s: distribution of strongly connected components
%   v: singular values
%   V: left and right singular vector
% Default value is '-p:cdhwsvV' which does all the plots
% 

% Transform the matrix into the Snap input graph format
tempSnapGraphfile = '../../Snap/temp.snap';                      %********************
saveMatrix2FromToGraphFile(graphMatrix, tempSnapGraphfile); 

eignvTop = 3; 

% -------------------  get common snap features --------------------------
% run Snap on the graph 
curDir = pwd;
cd '../../Snap';   %********************
cmd = '!SnapApp.exe'; 


% ---  for singular value  --
command = sprintf('%s -i:temp.snap -d:F -o:OUT -t:\" %s \" -p:v > OUT.snap.log ', cmd, 'curSubG' ) ;
eval(command)
inputfile = sprintf('sval.OUT.tab'); 
[rank, value ] = textread(inputfile, '%d  %f' , 'headerlines', 4); 
provideLenght = length(value); 

eigenVfeatures = zeros(eignvTop, 1); 
if ( provideLenght < eignvTop)
    eigenVfeatures(1:provideLenght) = value(1:provideLenght);
else
    eigenVfeatures(1:eignvTop) = value(1:eignvTop);    
end 

cd(curDir); 



if (opFlag)

    
% curDir = pwd;
% cd '../Snap'; 
% cmd = '!SnapApp.exe'; 
% command = sprintf('%s -i:temp.snap -d:F -o:OUT -t:\" %s \" > OUT.snap.log ', cmd, graphtitle ) 
% eval(command)
% cd(curDir); 
% 
%     
%     
% % process the output properties files 
% prefix = {'cDD';'dd';'hop';'lsv';'rsv'; 'scc'; 'sval'; 'wcc'}
% filesize = size(prefix, 1); 
% path = '../Snap/'; 
% for i = 1: filesize
%     inputfile = sprintf('%s%s.OUT.plt' , path, prefix{i}); 
%     outputfile = sprintf('%s.%s.plt', outputPre,  prefix{i}); 
%     cmd = sprintf('!cat %s > %s ', inputfile, outputfile) 
%     eval(cmd)
% 
%     inputfile = sprintf('%s%s.OUT.png' , path, prefix{i}); 
%     outputfile = sprintf('%s.%s.png', outputPre,  prefix{i}); 
%     cmd = sprintf('!cat %s > %s ', inputfile, outputfile) 
%     eval(cmd)    
% 
%     inputfile = sprintf('%s%s.OUT.tab' , path, prefix{i}); 
%     outputfile = sprintf('%s.%s.tab', outputPre,  prefix{i}); 
%     cmd = sprintf('!cat %s > %s ', inputfile, outputfile)  
%     eval(cmd)    
% end 
% 
% 
% 
% % -------------------  need to analyze the connectComponent file  --------------------------
% file = sprintf('%s.scc.tab' , outputPre ); 
% [ r, value ] = textread(file, '%d  %d', 'headerlines', 4); 
% largestccnodeSize = r(1)
% 
% 
% 
% 
% % -------------------  need to get hop plots multiple times to become reliable estimation --------------------------
% hopTimes = 10; 
% 
% % -p:h for only the hop plots 
% path = './Snap/'; 
% prefix = 'hop'; 
% effectivediameterArray = zeros( hopTimes, 1 ); 
% 
% hopfig = sprintf('%s.hopplots.fig', outputPre); 
% H = figure; 
% title(sprintf('%s hop plots (jump vs. pairs.visted) - repeat %d times ', graphtitle, hopTimes)); 
% xlabel(' hop Distance ')
% ylabel(' Number of pairs of nodes ')
% hold on 
% 
% for j = 1: hopTimes
%     cd './Snap'; 
%     cmd = '!SnapApp.exe'; 
%     command = sprintf('%s -i:temp.snap -d:F -p:h -o:OUT -t:\" %s \" > OUT.snap.log ', cmd, graphtitle ) 
%     eval(command)
%     cd ..; 
% 
%     inputfile = sprintf('%s%s.OUT.plt' , path, prefix ); 
%     outputfile = sprintf('%s.%s.%d.plt', outputPre,  prefix, j ); 
%     cmd = sprintf('!cat %s > %s ', inputfile, outputfile) 
%     eval(cmd)
% 
%     inputfile = sprintf('%s%s.OUT.png' , path, prefix); 
%     outputfile = sprintf('%s.%s.%d.png', outputPre,  prefix, j ); 
%     cmd = sprintf('!cat %s > %s ', inputfile, outputfile) 
%     eval(cmd)    
% 
%     inputfile = sprintf('%s%s.OUT.tab' , path, prefix ); 
%     outputfile = sprintf('%s.%s.%d.tab', outputPre,  prefix, j ); 
%     cmd = sprintf('!cat %s > %s ', inputfile, outputfile)  
%     eval(cmd)    
% 
%     [hop, cc ] = textread(inputfile, '%d  %f' , 'headerlines', 4); 
%     hold on; 
%     loglog(hop, cc, '+-');     
%     effectivediameterArray(j) = size(hop, 1);
% end 
% 
% effectivediameterArray'
% 
% saveas(H, hopfig, 'fig'); 

end 