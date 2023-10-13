%function HCP_SVR_analysis_nodes_null(beh_type,variant_form,variant_type,beh_var)
%function HCP_SVR_analysis_nodes()
%
% based on HCP_SVR_analysis_cg.m
%
% Components of script to allow for node-based version of function
% 
% CG 9.30.2022
% Heavily based on script/example from A. Nielsen
%
% beh_type = 'YeoBeh' or 'SeitzmanBeh'
% variant_form = 'network', 'location', 'locationXnetwork'
% variant_type = 'border','ectopic'

% beh_var = index for which of the behavioral variables to run
% 
% NOTE: this version does NULL permutations 

% set up variables from batch call
beh_type = input_arg1;
variant_form = input_arg2;
variant_type = input_arg3;
beh_var = str2double(input_arg4);

num_perms = 1000; % set not too high for now
beh_regress = true; % set up to regress out nuisance vars from data

% for debugging
% save('temp0_out.mat');

% set paths
top_dir = '/home/cgv5452/Ashley-behprediction/';
input_dir = [top_dir 'input_data_noBlueBrains/'];
output_dir = [top_dir 'output_data_noBlueBrains/' beh_type '_' variant_form '_prediction_results_byvar/'];
if beh_regress
    output_dir = [top_dir 'output_data_noBlueBrains_behregress/' beh_type '_' variant_form '_prediction_results_byvar/'];
end

% load subject data
subInfoTable = load_subject_data(input_dir);

% Load fMRI data 
thisVariantInput = load_fMRI_data(input_dir, variant_form, variant_type);

%debug
%save('temp1_out.mat');

% Load behavioral data
[behScores, behLabels] = load_behavioral_data(input_dir,beh_type);

if beh_regress    
    % regress out nuisance variables as in Kong...Yeo
    % vars: age, sex, FD, DVARS, BMI, and total brain volume (don't have
    % brain volume)
    % NOTE: should move this within CV since they "estimate regression
    % coefficients on training folds only and apply to test folds"    
    nuis_vars = [subInfoTable.age subInfoTable.sex subInfoTable.BMI subInfoTable.FD subInfoTable.DVARS];
    %[betas,temp2,resids] = regress(behScores',nuis_vars);
    %behScores = resids;
end

% remove subjects missing behavioral variables to make things cleaner
bhmean = mean(behScores,2); % take mean across behavioral variables to expose any nan
goodsubs = ~isnan(bhmean);
disp(['Good subs: ' num2str(sum(goodsubs))]);

% load pre-specified cross_validation groups (created through
% cvsets_nestfamily_updated.m)
load([input_dir '/HCP_CVsets_nestedfamily.mat']);

% SELECT VALUES USED (restrict to good subs and beahvioral variables)
% restrict only to subs with scores on all measures for ease
theseScores_init = behScores(goodsubs,beh_var);
thisVariantInput = thisVariantInput(:,:,goodsubs);
thisCVInput = subs_byCVnum(goodsubs);


% CROSSVALIDATION: 10-fold

% initialize variables
predict_r = nan; %*ones(length(behLabels),1);
predict_r_perfold = nan*ones(1,10);
ectopic_r_perfold = nan*ones(1,10);

% for checks
score_num_vals = numel(unique(theseScores_init));
num_nan_subs = sum(goodsubs); % should always be the same now
if score_num_vals > 2
    
    % debug
    %save('temp2_out.mat');
    
    rng(beh_var); % seed predictably for replicability
    
    for n = 1:num_perms
        
        disp(['XXXX NULL PERMUTATION #: ' num2str(n)]);
        
        %%%% FOR NULL PERMUTATIONS: PERMUTE SCORE ORDER
        theseScores = theseScores_init(randperm(length(theseScores_init)));
        
        if beh_regress
            %predict_info = tenFOLD_svm_scripts_matlab(thisVariantInput,theseScores,0,0,'SVR',thisCVInput,nuis_vars);
            %predict_info = tenFOLD_svm_scripts_matlab_v2(thisVariantInput,theseScores,0,0,'SVR',thisCVInput,nuis_vars);
            predict_info = tenFOLD_svm_Dworetsky(thisVariantInput,theseScores,0,0,'SVR',thisCVInput,nuis_vars);
        else
            %predict_info = tenFOLD_svm_scripts_matlab(thisVariantInput,theseScores,0,0,'SVR',thisCVInput);
            predict_info = tenFOLD_svm_Dworetsky(thisVariantInput,theseScores,0,0,'SVR',thisCVInput);
        end

        % store this info in an easier way for plotting       
        predict_r(n) = corr(theseScores,predict_info.predictedLabels);
        predict_r_perfold(n,:) = predict_info.rvalues_perfold;
        
        %debug
        %save('temp3_out.mat');
     
        clear predict_info;
    end
    
    % save out detailed results in case needed
    fout = [output_dir 'SVR_'  beh_type '_var' num2str(beh_var) '_' variant_form '_' variant_type '_prediction_results_NULL.mat'];
    save(fout,'predict_r','predict_r_perfold','behLabels','beh_var','num_nan_subs','score_num_vals');
    %fout_extended = [output_dir 'SVR_'  beh_type '_var' num2str(beh_var) '_' variant_form '_' variant_type '_prediction_results_extended_NULL.mat'];
    %save(fout_extended,'predict_info','-v7.3');
    
else
    disp(['WARNING: ' behLabels{beh_var} ' skipped due to too few possible values']);
    % Yeo 19 omitted because very little variance (almost all subs had the same
    % score)
    % Note that distribution on 39 and 58 also a bit weird, many not normal
end
    

