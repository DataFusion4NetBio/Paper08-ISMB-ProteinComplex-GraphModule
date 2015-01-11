%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
%%  copyright @  Yanjun Qi - qyj@cs.cmu.edu 
%%  
%%  Please cite: 
%%  Y. Qi, F. Balem, C. Faloutsos, J. Klein-Seetharaman, Z. Bar-Joseph, Protein Complex Identification by Supervised Graph Clustering,
%%  Bioinformatics 2008 24(13), i250-i268 (The 16th Annual International Conference Intelligent Systems for Molecular Biology (ISMB), July 2008
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% - trainMethodchoice {'svm', 'bc'}

function [] = complexTrainModel( trainCSVfile, derivedModelF, featureChoiceFlag, trainMethodchoice, trainMethodpara )

if ( nargin <5 )
    trainMethodchoice = 'svm'; 
    trainMethodpara =  '-t 0 -e 0.0001';  
end 


trainCSVfile = '../../temp/temp.train.csv'; 
temp_train_file = '../../temp/temp.train.data'; 
garbage_file = sprintf('%s.train.log', derivedModelF);



% ---- combine the feature files into one train file and convert to the format we need -------
cmd = '!perl ../trainTest/combineTrainPNfeature2csv.pl '; 
command = sprintf('%s %s %s %s ', cmd, TrainPfeaFile, TrainNfeaFile, trainCSVfile) 
eval(command); 


% ---- traing to get the model -------

if strcmp(trainMethodchoice, 'svm')
    % convert the two sets into SVML format
    cv_cmd = sprintf('!perl ../trainTest/convert_file_forSvm.pl %s %s', trainCSVfile, temp_train_file); 
    eval(cv_cmd); 

    % set up commands 
    learner = './trainTest/svm_learn.exe';
    train_cmd = sprintf('!%s -# 1000 %s %s %s > %s.train.log', learner, trainMethodpara, temp_train_file, derivedModelF, derivedModelF)
    eval(train_cmd)

elseif strcmp(trainMethodchoice, 'bc')
    % 'bc' could only apply on the {'29feasvd', 'densityNode'}
    learner = 'perl ../trainTest/bcTrain_12ppiNetv8.pl ';
    train_cmd = sprintf('!%s %s %s %s %s 1 > %s.train.log', learner, trainMethodpara, trainCSVfile, featureChoiceFlag, derivedModelF , derivedModelF)
    eval(train_cmd)    
end 
