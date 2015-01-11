%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---- test to get the prediction score of one cluster (feature array input) -------
%
% function [] = complexTestOneCluster( trainedModelF, testFeaOneCluster, predictedOuputF, featureChoiceFlag, trainMethodchoice )
%
%  - trainMethodchoice {'svm', 'bc'}

function [] = complexTestOneCluster( trainedModelF, testFeaOneCluster, predictedOuputF, featureChoiceFlag, trainMethodchoice )


if ( nargin <5 )
    trainMethodchoice = 'svm'; 
end 


garbage_file = '../../temp/temp.test.log';
testCSVfile = '../../temp/temp.test.csv'; 
temp_test_file = '../../temp/temp.test';



% ---- combine the feature file into one csv format data file -------
vfeature = testFeaOneCluster ; 
fid = fopen( testCSVfile ,'w');
fprintf(fid,'%.6f', vfeature(1) );
for j =2: length( vfeature )
    fprintf(fid,',%.6f', vfeature(j) );    
end 
fprintf(fid,',1\n'); 
fclose(fid); 




% ---- testing by the trained model -------
if strcmp(trainMethodchoice, 'svm')
    % convert the sets into SVML format
    cv_cmd = sprintf('!perl ../trainTest/convert_file_forSvm.pl %s %s > %s', testCSVfile, temp_test_file, garbage_file); 
    eval(cv_cmd); 

    % set up commands 
    classifier = 'e:\qyj\research\learning_tools\svmlight\v60\svm_classify.exe';
    test_cmd = sprintf('!%s %s %s %s > %s', classifier, temp_test_file, trainedModelF, predictedOuputF, garbage_file) ;
    eval(test_cmd); 
    
elseif strcmp(trainMethodchoice, 'bc')
     % 'bc' could only apply on the {'29feasvd', 'densityNode'}    
    classifier = 'perl ../trainTest/bcTest_12ppiNetv8.pl  ';
    test_cmd = sprintf('!%s %s %s %s %s > %s', classifier,  trainedModelF, featureChoiceFlag, testCSVfile, predictedOuputF, garbage_file) ;
    eval(test_cmd);         
end 