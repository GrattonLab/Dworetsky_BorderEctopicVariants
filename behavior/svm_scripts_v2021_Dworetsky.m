function [ results_struct ] = svm_scripts_v2021_Dworetsky( rmat, labels, numFeatures, testRmat, testLabels, type_of_prediction, option)

%Last Updated: 01.22.21 Ashley Nielsen
% 06/13/2023 Caterina Gratton
%   Edited script to simplify and clean up for release
%
% 09/28/2022 Caterina Gratton
%   Changed to input type of prediction
%
% 02.12.19 Ashley Nielsen
%   Added accompanying tenfoldCV_wrapper.m last located at
%   /data/cn5/ashley/tenfoldCV_wrapper.m
%
% 05.08.17 Ashley Nielsen
%   Added multi-class SVM
%   Added a test individuals on multi-class SVM
%
% 03.28.17 Ashley Nielsen
%   Added an option to not do LOOCV (useful for train/test sets)
%
% 08.18.16 Ashley Nielsen
%   Added the ability to not do feature selection
%   Made output a structure
%
% 07.22.16 Ashley Nielsen
%   Fixed issue if the input isn't square in variablefeatureselector
%
% 06.14.16 Ashley Nielsen
%   Fixed up the R to Z norm so that don't do this if input > 1
%   Made outputs struct instead of many variables
%   Added a re-calculation of SVM/SVR using only consensus features
%
% 03.02.16 Ashley Nielsen
%   Cleaned up the outputs so they are more useful

% This script can run SVM or SVR depending on the inputs. This script can
% also test a new dataset on models trained from another dataset.

% INPUTS:
%   rmat - features x numSubjects matrix of correlations
%           *In the normal SVM/SVR condition, these are the data you wish to
%           classify/predict label
%           *In the test SVM/SVR condition, these are the data you are
%           using to classify/predict the label of your test data
%
%   labels - numSubjects list of data
%           *In the SVR condition, this should be the variables that you
%           want to predict (ex: age, meanFD, tic score...)
%
%   numFeatures - list of the number of features to test ex: [200, 400, 600]
%           *If you would like to do NO feature selection, make this 0 =
%           OUR CASE
%
%   testRmat - features x number of test subjects matrix of correlations
%
%   testLabels - number of test subjects list of labels
%
%   type_of_prediction - type of prediction to do
%           ('SVM','SVR','multi-class')
%           CG: CODE currently simplied to just run SVR case
%
%   nus_regress - either empty or array of nuisance variables to regress from labels 
%
%   option - 2 if you do not want to do LOOCV = OUR CASE (simplified --
%   other options previously offered)

% OUTPUTS:
%   predictedLabels - numSubjects x length(numFeatures) matrix of predicted
%           class/continuous labels
%           *In the SVM condition, this will be a list of 1s and 0s which
%           will tell you whether that subject was predicted correctly (1 =
%           correct, 0 = incorrect) at each numFeatures
%           *In the SVR condition, this will be a list of the SVR predicted
%           labels for each subject at each numFeatures (ex: predicted age)

%   featureList - max(numFeatures) x 2 x numSubjects x length(numFeatures)
%           list of the ROI pairs determined for the SVM/SVR prediction

%   featureWeights - max(numFeatures) x numSubjects x length(numFeatures)
%           weights of features used by the SVM and SVR algorithms

%   performanceMeasure - various measures of performance
%           *In the SVM condition, this will give you the hitRate (%
%           correct) for each numFeatures tested
%           *In the SVR condition, this will give you the R^2 or the
%           variance explained by the these features
%           *In the SVM test subjects condition, this will do something...
%           *In the SVR test subjects condition, this will give you the
%           predicted labels for the test subjects


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize Variables
consensusFeatureMat = zeros(size(rmat,1),size(rmat,2),length(numFeatures));

num_vals_data = numel(unique(labels)); % still check this to confirm no errors
disp(['# values labels take: ' num2str(num_vals_data) ]);

% set up to work with original code (a bit clunky)
switch type_of_prediction
    case 'SVR'
        svr_check = true;
        svm_check = false;
        svm_multi_check = false;
        if num_vals_data < 10
            disp('WARNING: few values for SVR!');
        end
    otherwise
        error('code simplified to only do SVR for release')
end


% set up some defaults
[numRows numCols numSubjects] = size(rmat);
noFeatSelect = true; % we're not feature selecting in this version of the code
numFeatures = length(nonzeros(sum(rmat,3))); %assumes a vector input of features
numSubjects = length(labels);

% set up empty variables for filling
predictedLabels = zeros(length(testLabels),length(numFeatures));
featureList = zeros(max(numFeatures),2,numSubjects,length(numFeatures));

%Re-organize labels as code expects, if needed
if size(labels,1) == 1
    labels = labels';
end
if size(testLabels,1) == 1
    testLabels = testLabels';
end

%Convert R Correlations to Z values (Normalization)
if max(abs(rmat(:)))<=1 %max(max(max(abs(rmat)))) <= 1 %WTF why is 1 not equal to 1
    if length(unique(rmat(:))) ~=2
        rmat = atanh(rmat); %convertR2Z(rmat);
        testRmat = atanh(testRmat); %convertR2Z(testRmat);
    else
        disp('Binary values so skipping r to z transform...')
    end
else
    disp('Values greater than 1 so skipping r to z transform...')
end
%rmat(rmat==Inf)=1;


%% Main Sequence
%% SVR with test subjects

%Single Model in the Training Set
disp('Create Single Model in Training Set')
[featureList_temp,features] = variablefeatureSelector(rmat,labels,numFeatures,svr_check,svm_multi_check,noFeatSelect,option); % select all features

for f = 1:length(numFeatures) % if there's more than one feature set (not the case here)
    
    %Extract features from training data
    nonZFeatures = length(nonzeros(featureList_temp(:,1,f)));
    %%CG: why are we looking for non-zeros here? should always
    %%be non-zero since just the edge path NM I see - if we had
    %%selected above some might be non-zero
    featureList(1:nonZFeatures,:,f) = featureList_temp(1:nonZFeatures,:,f);

    trainFeatures = features(:,1:nonZFeatures,f);

    % fit the SVM model to training data
    mdl = fitrsvm(trainFeatures,labels);

    % set up the test features (clunky? but allows for more flexible
    % feature selection cases than what we use here)
    testsFeatures = zeros(length(testLabels),nonZFeatures);
    for t = 1:length(testLabels)
        for j = 1:nonZFeatures
            testsFeatures(t,j) = testRmat(featureList(j,1,f),featureList(j,2,f),t);
        end
    end

    % predict values in held out test data
    predictedLabels = predict(mdl,testsFeatures);
    rvalue = corr(testLabels,predictedLabels);

    % RESULTS - store
    results_struct(f).predictedLabels = predictedLabels;
    results_struct(f).featureList = featureList(:,:,:,f);
    %results_struct(f).featureWeights = featureWeights(:,:,f);
    results_struct(f).rvalue = rvalue;
    results_struct(f).R2 = rvalue.^2;

end



end

function [featureList, features] = variablefeatureSelector(rmat,labels,numFeatures,svr_check,svm_multi_check,noFeatSelect,option)
% previous more general version could do feature selection
% here, simplified to just piece that checks that features are > 0 in at
% least some people (comes up in location based analysis)

[numRows,numCols,numSubjects] = size(rmat);
testF = length(numFeatures);
if noFeatSelect % if you want to bypass feature selection
    nonZFeatureMat = sum(rmat,3)~=0;

    count = 1;
    for x = 1:numRows
        for y = 1:numCols
            if nonZFeatureMat(x,y) ~= 0 % CG added so zeros would be skipped
                featureList(count,:) = [x y];
                count = count+1;
            end
        end
    end
    for n = 1:length(featureList)
        features(:,n) = squeeze(rmat(featureList(n,1),featureList(n,2),:));
    end

end
end



