function results = tenFOLD_svm_Dworetsky(corrmat, labels, ind_test_set_corrmat, ind_test_set_labels,type_of_prediction,CVinput,varargin)
% function results = tenFOLD_svm_scripts_matlab(corrmat, labels, ind_test_set_corrmat, ind_test_set_labels,type_of_prediction,varargin)
%
% INPUTS

% corrmat - numROIs X numROIs X numSubjects
% labels - what you are trying to predict with FC - numSubjects X 1 - Ok if you included NaNs
        % could be BINARY: e.g., insurance status (0,1)
        % could be CATEGORICAL (up to 20): e.g., education (1,2,3,4)
        % could be CONTINUOUS: e.g., age (36.4-45.2 weeks)

% ind_test_set_corrmat - numROIs x numROIs X numTESTSubjects or 0 if you aren't going to test an independent set
% ind_test_set_labels - numTESTSubjects X 1 or 0 if you aren't going to test an independent  set
% type_of_prediction: 'SVR', 'SVM', 'multiclass'
% cv_sets - [OPTIONAL] cell array as long as numbr of CV runs (assume 10 for now),
%   each containing the test set IDS of subjects (if not provided, divide randomly)
%   note, IDs expressed as ID of corrmat structure
% nuis_vars - [OPTIONAL] array with nuisance variables to regress from labels as part of CV as in Kong et al. [CG]

% OUTPUTS - STRUCTURE w/ these variables

% predictedLabels - numSubjects X 1
        % for the SVM case: 1 if correct, 0 if incorrect
        % for the SVR case: predicted continuous label
        % for the MC-SVM case: class # a subject was categorized as (e.g., 3 out 1,2,3,4,5)
% R2 - 1 value; only for the SVR case, correlation between the true/predicted labels ^2
% hitRate - 1 value
        % for the SVM case: % of subjects correctly classified (if equal # in each class,  chance = 50%)
        % for the MC-SVM: % of subjects correctly classified (chance is < 50% and depends on #/size of classes)
% featureWeights
        % for the SVM/SVR case: {10 x 1} cell with numROIs x numROIS
            % each cell is the featureWeights from 1 fold of cross-validation
            % each matrix is how all the included FC was weighted in training 
                % featurWeights can be visualized w/ BLAH file I made for Jeanette
        % for the MC-SVM case: {10 x 1} cell with another layer of cells
            % feature weights for every possible combo of classes (e.g., 1vs2, 1vs3, 1vs4, 2vs3, 2vs4, 3vs4) 
% N - 1 value; number of subjects w/o NaNs
% predictedTestLabels - numTESTsubjects x 1
        % for the SVM case: 1 if correct, 0 if incorrect
        % for the SVR case: predicted continuous label
        % for the MC-SVM case: class # a subject was categorized as (e.g., 3 out 1,2,3,4,5)
% testR2 - 1 value; only for the SVR case, correlation between the true/predicted labels ^2 in the independent dataset
% testHitRate - 1 value for the independent datset
        % for the SVM case: % of subjects correctly classified (if equal # in each class,  chance = 50%)
        % for the MC-SVM: % of subjects correctly classified (chance is < 50% and depends on #/size of classes)


% DEPENDENCIES:
% svm_scripts_v2021_matlab.m

% USAGE:
% eLABE_predictAGE = tenFOLD_svm_scripts_matlab(corrmat,PMA,0,0);  % SVR on a single dataset
% predictAGE_trainELABE_testCUDDEL = tenFOLD_svm_scripts_matlan(corrmat_eLABE,PMA_eLABE,corrmat_CUDDEL,PMA_CUDDEL); % SVR w/ a independent test set
% eLABE_predictINSUR = tenFOLD_svm_scripts_matlab(corrmat,insur,0,0); % SVM on a single dataset
% eLABE_predictEDU = tenFOLD_svm_scripts_matlab(corrmat,edu,0,0); % MC-SVM on a single datasetn
%
% Original script by A. Nielsen
% 09.27.2022 CG - made a few clarification edits/saving edits
%                 changed to make SVR/SVM/multiclass an input and to pass
%                 it along
% 06.13.2023 CG - cleaned up script options substantially for readability


% IDENTIFY ONLY SUBJECTS w/ no NaNs
idx_nonNaN = find(~isnan(labels)==1);
numSubjects = length(labels(idx_nonNaN));

% CREATE THE 10-FOLD TRAINING and TEST SETS
% cg - edited to feed in CV division (pre-specified to keep familial
% structure in account)
% input CV set division (first varargin)
for n = 1:1
    for cv = 1:10
        test_set{n,cv} = find(CVinput==cv); %for now only have one cv set
        train_set{n,cv} = setdiff(1:numSubjects,test_set{n,cv});
    end
end

% look for varargin with nuisance variable information to pass on
if length(varargin) > 0
    nuis_vars = varargin{1};
else
    nuis_vars = [];
end


numIterations = size(train_set,1);
if size(labels,1)==1
    labels = labels';
end

if size(ind_test_set_labels,1)==1
    ind_test_set_labels = ind_test_set_labels';
end




% WHAT TYPE OF SVM ARE WE DOING HERE?
% CG - changed this to an input
switch type_of_prediction
    case 'SVR'
        cont_svr = true;
        binary_svm = false;
        multi_class = false;
    case 'SVM'
        cont_svr = false;
        binary_svm = true;
        multi_class = false;
    case 'multi_class'
        cont_svr = false;
        binary_svm = false;
        multi_class = true;
end
        

% ENTER INTO FOLDS OF TRAINING/TESTING
for n = 1:numIterations
    disp(['Iter #',num2str(n)])
    % initialize a few variables
    predictedLabels = nan.*ones(length(labels),1);
    predictedTestLabels = nan.*ones(length(ind_test_set_labels),10);
    hitRate = zeros(10,1);
    
    for cv = 1:10 % loop throug cross-validation folds
        disp(['FOLD #',num2str(cv)])
        
        % adjust labels to account for nuisance regression
        train_labels = labels(idx_nonNaN(train_set{n,cv}));
        test_labels = labels(idx_nonNaN(test_set{n,cv}));
        if ~isempty(nuis_vars)
            %[betas_nuis_train,temp2,resids] = (train_labels,nuis_vars(train_set{n,cv}));
            mdl_nuis_train = fitlm(nuis_vars(train_set{n,cv},:),train_labels);
            train_labels_fin = train_labels - predict(mdl_nuis_train,nuis_vars(train_set{n,cv},:));
            test_labels_fin = test_labels - predict(mdl_nuis_train,nuis_vars(test_set{n,cv},:));
        else
            train_labels_fin = train_labels;
            test_labels_fin = test_labels;
        end
        
        % input data: 
        trainFeatures = corrmat(:,:,idx_nonNaN(train_set{n,cv}));
        testFeatures = corrmat(:,:,idx_nonNaN(test_set{n,cv}));
        
        % TRAINING/TESTING
        temp = svm_scripts_v2021_Dworetsky(trainFeatures,train_labels_fin,0,testFeatures,...
            test_labels_fin,type_of_prediction,2);
        
        % REORGANIZE PREDICTED LABELS FROM TEST SET BACK into the TOTAL SET
        for i = 1:length(idx_nonNaN(test_set{n,cv})) % For SVM, SVR, and MULTI-CLASS
            predictedLabels(idx_nonNaN(test_set{n,cv}(i))) = temp.predictedLabels(i);
        end
        
        % hold onto R^2 info for my perusal - CG
        r2_eachfold(n,cv) = temp.R2;
        r_eachfold(n,cv) = temp.rvalue;
        
        % INDEPENDENT TEST SET
        if size(ind_test_set_corrmat,1) > 1 
            temp = svm_scripts_v2021_matlab(corrmat(:,:,idx_nonNaN(train_set{n,cv})),labels(idx_nonNaN(train_set{n,cv})),0,ind_test_set_corrmat,ind_test_set_labels,2);
            
            predictedTestLabels(:,cv) = temp.predictedLabels;   % For SVM, SVR, and MULTI-CLASS
            
            if binary_svm || multi_class
                testHitRate(cv) = temp.hitRate;
            end
        end
        
%         % SAVE OUT THE FEATURE WEIGHTS PER FOLD
%         if ~multi_class
%             fW = zeros(size(corrmat,1));
%             for i = 1:size(temp.featureList,1)
%                 fW(temp.featureList(i,1,1),temp.featureList(i,2,1)) = temp.featureWeights(i,1);
%                 fW(temp.featureList(i,2,1),temp.featureList(i,1,1)) = temp.featureWeights(i,1);
%             end
%             featureWeights{cv} = fW;
%         else
%             featureWeights{cv} = temp.featureWeights;
%         end
%     end
%     

    end   % SAVE RESULTS INTO A STRUCTURE

    results(n).predictedLabels = predictedLabels;
    %     results(n).featureWeights = featureWeights;

    if cont_svr
        results(n).R2 = corr(labels(idx_nonNaN),predictedLabels(idx_nonNaN)).^2;
        results(n).N = length(idx_nonNaN);
        results(n).R2_perfold = r2_eachfold(n,:);
        results(n).rvalues_perfold = r_eachfold(n,:);
    end
    if multi_class
        results(n).hitRate = mean(predictedLabels(idx_nonNaN)==labels(idx_nonNaN));
        results(n).N = length(idx_nonNaN);
    end
    if binary_svm
        results(n).hitRate = mean(predictedLabels(idx_nonNaN));
        results(n).N = length(idx_nonNaN);
    end

    % SAVE INDEPENDENT TEST SET RESULTS INTO A STRUCTURE
    if size(ind_test_set_corrmat,1)>1
        results(n).predictedTestLabels = predictedTestLabels;

        if cont_svr
            results(n).testR2 = (corr(ind_test_set_labels,mean(predictedTestLabels,2))).^2;
        end
        if binary_svm
            results(n).testHitRate = mean(testHitRate);
        end
        if multi_class
            results(n).testHitRate = mean(predictedTestLabels==ind_test_set_labels);
        end
    end



end