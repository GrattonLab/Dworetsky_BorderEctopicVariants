%function HCP_SVR_nodes_finalplotting()
% function HCP_SVR_nodes_finalplotting()
% 
% Plot out svr data calculated through nodes loop

% beh_type = 'YeoBeh' or 'SeitzmanBeh'
% variant_form = 'network', 'location', 'locationXnetwork'
% variant_type = 'border','ectopic'
% beh_var = index for which of the behavioral variables to run


clear all; close all;

%%% set up variables 
beh_type = 'YeoBeh';
variant_form = 'location';
compare_to_null = false; % set to false if not run yet
beh_regress = true; % were other variables regressed from behavior
parc_free = true; %run on the parcellation free variants or originals?

%%% set paths
top_dir = '/home/cgv5452/Ashley-behprediction/';
if ~parc_free
    input_dir = [top_dir 'input_data_noBlueBrains/'];
    output_dir = [top_dir 'output_data_noBlueBrains/' beh_type '_' variant_form '_prediction_results_byvar/'];
    if beh_regress
        output_dir = [top_dir 'output_data_noBlueBrains_behregress/' beh_type '_' variant_form '_prediction_results_byvar/'];
    end
else
    % Parcellation Free version
    input_dir = [top_dir 'input_data_noBlueBrains_parcellationfree/'];
    output_dir = [top_dir 'output_data_noBlueBrains_parcellationfree/' beh_type '_' variant_form '_prediction_results_byvar/'];
    if beh_regress
        output_dir = [top_dir 'output_data_noBlueBrains_behregress_parcellationfree/' beh_type '_' variant_form '_prediction_results_byvar/'];
    end
end

%%% code
% load behavioral info
[behScores, behLabels] = load_behavioral_data(input_dir,beh_type);
behCategories = readtable('BehCategories.txt');
%behCategories = readtable('BehCategories2.txt');
behCatData = behCategories.Var1;
behCatNames = {'cognitive', 'soc&affective', 'lifestyle', 'sensory-motor', 'personality','other'};
%behCatNames = {'NIH-cognition','NIH-emotion','beh_test','MRI_test','self-report','NIH-sensmot'};

% remove subjects missing behavioral variables to make things cleaner
bhmean = mean(behScores,2); % take mean across behavioral variables to expose any nan
goodsubs = ~isnan(bhmean);
disp(['Good subs: ' num2str(sum(goodsubs))]);

% load data calculated through nodes script
border_r = nan*ones(length(behLabels),1);
ectopic_r = nan*ones(length(behLabels),1);
both_r = nan*ones(length(behLabels),1);
border_r_perfold = nan*ones(length(behLabels),10);
ectopic_r_perfold = nan*ones(length(behLabels),10);
both_r_perfold = nan*ones(length(behLabels),10);
score_num_vals = nan*ones(length(behLabels),1);
num_nan_subs = nan*ones(length(behLabels),1);
for b = 1:length(behLabels)
    
    % should have variables: 'predict_r','predict_r_perfold','behLabels','beh_var','num_nan_subs','score_num_vals'
    border_data_file = [output_dir 'SVR_'  beh_type '_var' num2str(b) '_' variant_form '_border_prediction_results.mat'];
    if exist(border_data_file)
        border_data = load(border_data_file);
        border_r(b) = border_data.predict_r;
        border_r_perfold(b,:) = border_data.predict_r_perfold;
        
        ectopic_data_file = [output_dir 'SVR_'  beh_type '_var' num2str(b) '_' variant_form '_ectopic_prediction_results.mat'];
        ectopic_data = load(ectopic_data_file);
        ectopic_r(b) = ectopic_data.predict_r;
        ectopic_r_perfold(b,:) = ectopic_data.predict_r_perfold;
        
        % control: SVR run on both together
        both_data_file = [output_dir 'SVR_'  beh_type '_var' num2str(b) '_' variant_form '_all_prediction_results.mat'];
        both_data = load(both_data_file);
        both_r(b) = both_data.predict_r;
        both_r_perfold(b,:) = both_data.predict_r_perfold;
        
        % for checks
        score_num_vals(b) = border_data.score_num_vals; %should be the same for border and ectopic
        num_nan_subs(b) = border_data.num_nan_subs; %should be the same for border and ectopic
        
        clear border_data; clear ectopic_data;
        
        if compare_to_null
            NULL_border_data_file = [output_dir 'SVR_'  beh_type '_var' num2str(b) '_' variant_form '_border_prediction_results_NULL.mat'];
            %NULL_border_data_file = [output_dir 'SVR_'  beh_type '_var1_' variant_form '_ectopic_prediction_results_NULL.mat'];
            NULL_border_data = load(NULL_border_data_file);
            NULL_border_r(b,:) = NULL_border_data.predict_r;
            NULL_border_r_perfold(b,:,:) = NULL_border_data.predict_r_perfold;
            
            NULL_ectopic_data_file = [output_dir 'SVR_'  beh_type '_var' num2str(b) '_' variant_form '_ectopic_prediction_results_NULL.mat'];
            %NULL_ectopic_data_file = [output_dir 'SVR_'  beh_type '_var1_' variant_form '_ectopic_prediction_results_NULL.mat'];
            NULL_ectopic_data = load(NULL_ectopic_data_file);
            NULL_ectopic_r(b,:) = NULL_ectopic_data.predict_r;
            NULL_ectopic_r_perfold(b,:,:) = NULL_ectopic_data.predict_r_perfold;
        end
        
        
    else
        disp(['WARNING: Var ' num2str(b) ' file does not exist... skipping']);
    end
    
end

% for comparisons, avg per fold (fisher transform first for normality, then
% undo)
avg_border_r_perfold = tanh(mean(atanh(border_r_perfold),2));
avg_ectopic_r_perfold = tanh(mean(atanh(ectopic_r_perfold),2));
avg_both_r_perfold = tanh(mean(atanh(both_r_perfold),2));

% COMPARE TO null
if compare_to_null
    nperms = size(NULL_border_r_perfold,2);
    avg_NULL_border_r_perfold = tanh(mean(atanh(NULL_border_r_perfold),3));
    avg_NULL_ectopic_r_perfold = tanh(mean(atanh(NULL_ectopic_r_perfold),3));
    per95_NULL_border = prctile(avg_NULL_border_r_perfold',95);
    per95_NULL_ectopic = prctile(avg_NULL_ectopic_r_perfold',95);

    % calculate whether border on AVG > null, etc.
    avgtot_border_r = tanh(nanmean(atanh(avg_border_r_perfold)));
    avgtot_NULLborder_r = tanh(nanmean(atanh(avg_NULL_border_r_perfold),1));
    pavg_border = (sum(avgtot_NULLborder_r > avgtot_border_r)+1)/(nperms+1); %add 1 to num and denom to make unbiased (?)
    
    avgtot_ectopic_r = tanh(nanmean(atanh(avg_ectopic_r_perfold)));
    avgtot_NULLectopic_r = tanh(nanmean(atanh(avg_NULL_ectopic_r_perfold),1));
    pavg_ectopic = (sum(avgtot_NULLectopic_r > avgtot_ectopic_r)+1)/(nperms+1); %add 1 to num and denom to make unbiased (?)
    disp(sprintf('Border mean r=%0.04f, p<%0.04f',avgtot_border_r,pavg_border));
    disp(sprintf('Ectopic mean r=%0.04f, p<%0.04f',avgtot_ectopic_r,pavg_ectopic));
    
    % make a figure of these results
    figure; 
    scatter(ones(length(avgtot_NULLborder_r),1),avgtot_NULLborder_r,'k.','jitter','on');
    hold on;
    plot(1,avgtot_border_r,'r.','MarkerSize',20);
    scatter(2*ones(length(avgtot_NULLectopic_r),1),avgtot_NULLectopic_r,'k.','jitter','on');
    plot(2,avgtot_ectopic_r,'r.','MarkerSize',20);
    xlim([0 3]);
    set(gca,'XTick',[1,2],...
      'XTickLabels',{'border','ectopic'},...
        'FontSize',8,...
        'FontName','Times');
    hline_new(0,'k-',1);
    ylabel('Average prediction (r)');
    title(['Avg behavioral prediction from variant ' variant_form ' vs. null']);
    save_fig(gcf,[output_dir 'SVR_' variant_form '_' beh_type '_avgpredict_vs_null.png']);

    % calculate p-values against the null for individual measures
    for b = 1:length(behLabels)
        pval_border(b) = (sum(avg_NULL_border_r_perfold(b,:) > avg_border_r_perfold(b))+1)/(nperms+1); %NOTE: need more perms to deal with MC 
        pval_ectopic(b) = (sum(avg_NULL_ectopic_r_perfold(b,:) > avg_ectopic_r_perfold(b))+1)/(nperms+1);
    end
    % FDR correction for p-values ... not sure if we want this since
    % correlated values, but for quick check for individual variables.
    % RERUN after more perms - p values not fine enough for now for this
    [hborder, tmp, tmp2, pFDR_border] = fdr_bh(pval_border,.05);
    [hectopic, tmp, tmp2, pFDR_ectopic] = fdr_bh(pval_ectopic,.05);    
end

% store information in a more digestable format
resultsTable = table(behLabels',num_nan_subs,score_num_vals,avg_border_r_perfold,avg_ectopic_r_perfold,'VariableNames',{'BehVar','N','numElem','BorderR','EctopicR'});
writetable(resultsTable,[output_dir 'SVR_' variant_form '_' beh_type '_prediction_results_avgrperfold.csv']);

% plot by category
cat_IDs = [1,2,4,5]; %categories to plot; none exist in lifestyle or other
figure; 
set(gcf,'position',[10 10 800 600]);
hold on;
for i = 1:length(cat_IDs)
    if compare_to_null
        scatter(i*ones(size(avg_NULL_border_r_perfold,2),1),nanmean(avg_NULL_border_r_perfold(behCatData == cat_IDs(i),:),1),'k.','jitter','on');
    end
    plot(i,nanmean(avg_border_r_perfold(behCatData == cat_IDs(i))),'b.','MarkerSize',30);
    plot(i,nanmean(avg_ectopic_r_perfold(behCatData == cat_IDs(i))),'r.','MarkerSize',30);
    
    if compare_to_null
        pval_category_border(i) = (sum(nanmean(avg_NULL_border_r_perfold(behCatData == cat_IDs(i),:),1) > nanmean(avg_border_r_perfold(behCatData == cat_IDs(i))))+1)/(nperms+1); %NOTE: need more perms to deal with MC 
        pval_category_ectopic(i) = (sum(nanmean(avg_NULL_ectopic_r_perfold(behCatData == cat_IDs(i),:),1) > nanmean(avg_ectopic_r_perfold(behCatData == cat_IDs(i))))+1)/(nperms+1); %NOTE: need more perms to deal with MC 
    end
end
xlim([0 length(cat_IDs)+1]);
set(gca,'XTick',[1:length(cat_IDs)],...
    'XTickLabels',behCatNames(cat_IDs),...
    'FontSize',9,...
    'FontName','Times');
ylabel('average prediction (r)');
save_fig(gcf,[output_dir 'SVR_' variant_form '_' beh_type '_category_predictionplot.png']);

% plot linearly, maybe against the null
[svals, sortedOrder] = sort(avg_ectopic_r_perfold,'ascend');
if strcmp(beh_type,'YeoBeh') && strcmp(variant_form,'location')
    [svals, sortedOrder] = sort(avg_border_r_perfold,'ascend');
    sortedOrder = sortedOrder(1:end); %remove nan
end
figure();
set(gcf,'position',[10 10 800 20*(length(sortedOrder))+200]);
plot(avg_border_r_perfold(sortedOrder),'b.','MarkerSize',10);
hold on;
plot(avg_ectopic_r_perfold(sortedOrder),'r.','MarkerSize',10);
plot(avg_both_r_perfold(sortedOrder),'g.','MarkerSize',10);
if compare_to_null
    plot(per95_NULL_border(sortedOrder),'--','Color',[.5 .5 1]);
    plot(per95_NULL_ectopic(sortedOrder),'--','Color',[1 .5 .5]);
%     for b = 1:size(avg_NULL_border_r_perfold,1)
%         swarmchart(b*ones(1,size(avg_NULL_border_r_perfold,2)),avg_NULL_border_r_perfold(b,:),5,[1 1 .9]);
%     end
end
legend({'border','ectopic','both','b-perm 95%','e-perm 95%'},'location','northwest');
ylabel('Prediction accuracy (r)');
title(['Prediction by variant: ' variant_form]);
xlim([0 length(behLabels)+1]);
set(gca,'XTick',[1:length(sortedOrder)],...
    'XTickLabels',behLabels(sortedOrder),...
    'FontSize',8,...
    'FontName','Times');
yline(0,'k-')
%hline_new(0.08,'k--',1)
view([90 -90]);
save_fig(gcf,[output_dir 'SVR_' variant_form '_' beh_type '_allpredictionplot.png']);



% Write some code for plotting results
figure;
scatter(border_r,ectopic_r);
xlabel('border prediction');
ylabel('ectopic prediction');
%hline_new(0,'k',1)
%vline_new(0,'k',1)
yline(0,'k');
xline(0,'k');
corr(border_r(~isnan(border_r)),ectopic_r(~isnan(ectopic_r)))
title('variant prediction levels');
save_fig(gcf,[output_dir 'SVR_' variant_form '_' beh_type '_borderXectopic_prediction_rtot.png']);

figure;
plot(avg_border_r_perfold,avg_ectopic_r_perfold,'k.','MarkerSize',20);
plot_min = -0.15; plot_max = 0.2;
xlim([plot_min,plot_max]);
ylim([plot_min,plot_max]);
axis square;
xlabel('border prediction');
ylabel('ectopic prediction');
if compare_to_null
    %patch([plot_min nanmean(per95_NULL_border) nanmean(per95_NULL_border) plot_min],[plot_min,plot_min,plot_max,plot_max],[0.8 .8 1],'FaceAlpha',0.25,'EdgeColor','none');
    %patch([plot_min plot_max plot_max plot_min],[plot_min,plot_min,nanmean(per95_NULL_ectopic),nanmean(per95_NULL_ectopic)],[1 .8 0.8],'FaceAlpha',0.25,'EdgeColor','none');
    vline_new(nanmean(per95_NULL_border),'b--',1);
    hline_new(nanmean(per95_NULL_ectopic),'r--',1);
end
%hline_new(0,'k',0.5);
%vline_new(0,'k',0.5);
yline(0,'k');
xline(0,'k');
corr(avg_border_r_perfold(~isnan(border_r)),avg_ectopic_r_perfold(~isnan(ectopic_r)))
title('variant prediction levels per fold');
save_fig(gcf,[output_dir 'SVR_' variant_form '_' beh_type '_borderXectopic_prediction_rfold.png']);

figure;
subplot(1,2,1)
plot(avg_border_r_perfold,avg_both_r_perfold,'k.','MarkerSize',20);
plot_min = -0.15; plot_max = 0.2;
xlim([plot_min,plot_max]);
ylim([plot_min,plot_max]);
hold on;
plot([plot_min plot_max],[plot_min, plot_max],'k');
axis square;
xlabel('border prediction');
ylabel('border+ectopic prediction');
yline(0,'k');
xline(0,'k');
corr(avg_border_r_perfold(~isnan(border_r)),avg_both_r_perfold(~isnan(ectopic_r)))
title('variant prediction levels per fold');

subplot(1,2,2)
plot(avg_ectopic_r_perfold,avg_both_r_perfold,'k.','MarkerSize',20);
plot_min = -0.15; plot_max = 0.2;
xlim([plot_min,plot_max]);
ylim([plot_min,plot_max]);
hold on;
plot([plot_min plot_max],[plot_min, plot_max],'k');
axis square;
xlabel('ectopic prediction');
ylabel('border+ectopic prediction');
yline(0,'k');
xline(0,'k');
corr(avg_ectopic_r_perfold(~isnan(border_r)),avg_both_r_perfold(~isnan(ectopic_r)))
title('variant prediction levels per fold');

save_fig(gcf,[output_dir 'SVR_' variant_form '_' beh_type '_borderXboth_prediction_rfold.png']);


figure;
subplot(1,2,1);
scatter(border_r,avg_border_r_perfold);
hold on;
plot([-0.2,0.2],[-0.2 0.2], 'k-');
%hline_new(0,'k',1);
%vline_new(0,'k',1);
yline(0,'k');
xline(0,'k');
axis square;
xlabel('r tot');
ylabel('r per fold');
title('border');

subplot(1,2,2);
scatter(ectopic_r,avg_ectopic_r_perfold);
hold on;
plot([-0.2,0.2],[-0.2 0.2], 'k-');
%hline_new(0,'k',1);
%vline_new(0,'k',1);
yline(0,'k');
xline(0,'k');
axis square;
xlabel('r tot');
ylabel('r per fold');
title('ectopic');
save_fig(gcf,[output_dir 'SVR_' variant_form '_' beh_type '_rtotXrfold.png']);
