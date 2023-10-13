function plot_task_activs_by_network(varargin)

% output directory is optional input - if not given, figures will appear but not save out


%% Setup info on subjects and contrasts
subs = dlmread('/gpfs/research/grattonlab/HCP/Tasks/analysis/HCPsubjects_for_contrasts.txt');
reassigned_dir = '/gpfs/research/grattonlab/HCP/Tasks/analysis/';
borderectopic_dir = '/gpfs/research/grattonlab/HCP/Tasks/analysis/';
contrastFilesLoc = '/gpfs/research/grattonlab/HCP/Tasks/analysis';
networks = {'DMN','Visual','FP','DAN','Sal','CO','SMd','SMl','Aud','PMN','PON'};
netIDs = [1 2 3 5 8 9 10 11 12 15 16]; %network IDs corresponding to 1-#networks
net_select = [1 2 3 4 6 9 11]; %selected networks (indices of netIDs) for analyses later
num_cortical_verts = 59412;

%% File info on variants
variantFile = [reassigned_dir '/HCPsubjects_task_reassigned.dscalar.nii'];
variant_ciftistruct = ft_read_cifti_mod(variantFile);
variantData = variant_ciftistruct.data;
border_variantData = variantData;
ectopic_variantData = variantData;
borderEctopicFile = [borderectopic_dir '/HCPsubjects_task_borderEctopic.dscalar.nii'];
borderectopic_ciftistruct = ft_read_cifti_mod(borderEctopicFile);
borderectopicData = borderectopic_ciftistruct.data;
border_variantData(borderectopicData == 2) = 0;
ectopic_variantData(borderectopicData == 1) = 0;


%% Set up contrast file data

%task_names = {'EMOTION' 'GAMBLING' 'LANGUAGE' 'MOTOR' 'RELATIONAL' 'SOCIAL' 'WM'};

contrasts = {'EMOTION_FACES-SHAPES' 'EMOTION_FACES' 'EMOTION_SHAPES' ...
    'GAMBLING_PUNISH-REWARD' 'GAMBLING_PUNISH' 'GAMBLING_REWARD' ...
    'LANGUAGE_MATH-STORY' 'LANGUAGE_MATH' 'LANGUAGE_STORY' ...
    'MOTOR_LF-AVG' 'MOTOR_LH-AVG' 'MOTOR_RF-AVG' 'MOTOR_RH-AVG' 'MOTOR_T-AVG' 'MOTOR_AVG' ...
    'RELATIONAL_MATCH-REL' 'RELATIONAL_MATCH' 'RELATIONAL_REL' ...
    'SOCIAL_RANDOM-TOM' 'SOCIAL_RANDOM' 'SOCIAL_TOM' ...
    'WM_2BK-0BK' 'WM_BODY-AVG' 'WM_FACE-AVG' 'WM_PLACE-AVG' 'WM_TOOL-AVG' 'WM_2BK' 'WM_0BK'};

contrasts_short = {'faces-shapes' 'faces' 'shapes' ...
    'punish-reward' 'punish' 'reward' ...
    'math-story' 'math' 'story' ...
    'lf-avg' 'lh-avg' 'rf-avg' 'rh-avg' 't-avg' 'avg' ...
    'match-rel' 'match' 'rel' ...
    'random-tom' 'random' 'tom' ...
    '2bk-0bk' 'body-avg' 'face-avg' 'place-avg' 'tool-avg' '2bk' '0bk'};

contrastFiles = cell(length(contrasts),1);
for c=1:length(contrasts)
    contrastFileName = sprintf('%s/HCPsubjects_%s.dscalar.nii', contrastFilesLoc, contrasts{c});
    temp = ft_read_cifti_mod(contrastFileName);
    contrastFiles{c} = temp.data;
    clear temp
end

%% File info on canonical network assignments
groupNetFile = '/gpfs/research/grattonlab/HCP/Tasks/analysis/HCP_infoMap_recolored.dtseries.nii';
groupNet_ciftistruct = ft_read_cifti_mod(groupNetFile);
groupNet = groupNet_ciftistruct.data(1:num_cortical_verts); %cortical data only

%% Do calculations

act_mean_border = nan(length(contrasts), length(subs), length(netIDs));
act_mean_ectopic = nan(length(contrasts), length(subs), length(netIDs));
act_canonical_mean = nan(length(contrasts), length(subs), length(netIDs));
act_canonical_other_mean = nan(length(contrasts), length(subs), length(netIDs));

for c = 1:length(contrasts)

    for s = 1:length(subs)
        contrastData = contrastFiles{c}(:, s); %already pre-selected cortex-only datapoints
        
        border_variantData_sub = border_variantData(:,s);
        ectopic_variantData_sub = ectopic_variantData(:,s);
        variantData_sub = variantData(:,s);

        for n = 1:length(networks)

            act_mean_border(c,s,n) = mean(contrastData(border_variantData_sub == netIDs(n)));
            act_mean_ectopic(c,s,n) = mean(contrastData(ectopic_variantData_sub == netIDs(n)));
            
            % do the same calculations for canonical regions in each network
            canonical_mask = groupNet == netIDs(n);
            canonical_mask(variantData_sub>0) = 0; % remove variants from canonical mask
            act_canonical_mean(c,s,n) = mean(contrastData(canonical_mask));
            
            % and do this for the average of canonical regions in all OTHER networks (i.e., all nonDMN regions)
            canonical_mask_other = groupNet ~= netIDs(n);
            canonical_mask_other(variantData_sub>0) = 0; % remove variants from this mask
            act_canonical_other_mean(c,s,n) = mean(contrastData(canonical_mask_other));
        end
    end
end

%% go by network by contrast and make plots
for n = 1:length(net_select)
    
    % set the minimum difference of abs(canonical - noncanonical)
    thresh = 0.5;

    % sort difference values in order of greatest to least (for plotting)
    net_ind= net_select(n);
    can_minus_noncan= zeros(length(contrasts),1);
    can_minus_noncan = nanmean(act_canonical_mean(:, :, net_ind), 2) - nanmean(act_canonical_other_mean(:, :, net_ind), 2);
    [sorted_diffs, sorti] = sort(can_minus_noncan, 'descend');
    sorted_diffs_highContrast = sorted_diffs; sorted_diffs_highContrast(abs(sorted_diffs) < thresh) = [];
    num_pos_contrasts = sum(sorted_diffs_highContrast > 0);
    sorti(abs(sorted_diffs) < thresh) = [];
    num_contrasts = length(sorti);
    contrast_names_short = contrasts_short(sorti);
    
    % find number of subjects with variant of that form and network
    % (just using first contrast to calculate, b/c same for all)
    num_border_subs = sum(~isnan(act_mean_border(1,:, net_ind)));
    num_ectopic_subs = sum(~isnan(act_mean_ectopic(1,:, net_ind)));
    
    % plot all 4 values for each contrast
    f=figure; hold on;
    errorbar((1:num_contrasts)-0.15,nanmean(act_mean_border(sorti, :, net_ind),2), ...
        nanstd(act_mean_border(sorti, :, net_ind),[],2)/sqrt(num_border_subs),'bo','LineWidth',2,'markersize',5); hold on;
    errorbar((1:num_contrasts)-0.05,nanmean(act_mean_ectopic(sorti, :, net_ind),2), ...
        nanstd(act_mean_ectopic(sorti,:,net_ind),[],2)/sqrt(num_ectopic_subs),'ro','LineWidth',2,'markersize',5); hold on;
    errorbar((1:num_contrasts)+0.05,nanmean(act_canonical_mean(sorti, :, net_ind),2), ...
        nanstd(act_canonical_mean(sorti,:,net_ind),[],2)/sqrt(length(subs)),'ko','LineWidth',2,'markersize',5); hold on;
    errorbar((1:num_contrasts)+0.15,nanmean(act_canonical_other_mean(sorti, :, net_ind),2), ...
        nanstd(act_canonical_other_mean(sorti,:,net_ind),[],2)/sqrt(length(subs)), 'o','color',[0.5 0.5 0.5],'LineWidth',2,'markersize',5); hold on;
    xlim([0.5 num_contrasts+0.5]);
    hline(0, 'k:');
    set(gca,'XTick',1:num_contrasts,'XTicklabel', contrast_names_short);
    set(gcf,'color','w');
    ylabel('avg. task activation');
    title(networks{net_ind},'fontsize',26);

    % change y-axis limits if max/min is higher/lower
    mins= [min(nanmean(act_mean_border(sorti, :, net_ind),2)), min(nanmean(act_mean_ectopic(sorti, :, net_ind),2)), ...
        min(nanmean(act_canonical_mean(sorti, :, net_ind),2)), min(nanmean(act_canonical_other_mean(sorti, :, net_ind),2))];
    minStds= [min(nanstd(act_mean_border(sorti, :, net_ind),[],2)/sqrt(length(subs))), ...
        min(nanstd(act_mean_ectopic(sorti, :, net_ind),[],2)/sqrt(length(subs))), ...
        min(nanstd(act_canonical_mean(sorti, :, net_ind),[],2)/sqrt(length(subs))), ...
        min(nanstd(act_canonical_other_mean(sorti, :, net_ind),[],2)/sqrt(length(subs)))];
    mins2 = mins - minStds; minVal = min(mins2);
    maxes= [max(nanmean(act_mean_border(sorti, :, net_ind),2)), max(nanmean(act_mean_ectopic(sorti, :, net_ind),2)), ...
        max(nanmean(act_canonical_mean(sorti, :, net_ind),2)), max(nanmean(act_canonical_other_mean(sorti, :, net_ind),2))];
    maxStds= [max(nanstd(act_mean_border(sorti, :, net_ind),[],2)/sqrt(length(subs))), ...
        max(nanstd(act_mean_ectopic(sorti, :, net_ind),[],2)/sqrt(length(subs))), ...
        max(nanstd(act_canonical_mean(sorti, :, net_ind),[],2)/sqrt(length(subs))), ...
        max(nanstd(act_canonical_other_mean(sorti, :, net_ind),[],2)/sqrt(length(subs)))];
    maxes2 = maxes + maxStds; maxVal = max(maxes2);
    if minVal<-1.5
        newMinVal = floor(minVal) + floor((minVal-floor(minVal))/0.05) * 0.05;
    else 
        newMinVal = -1.5;
    end
    if maxVal>3
        newMaxVal = floor(maxVal) + ceil((maxVal - floor(maxVal))/0.05) * 0.05;
    else
        newMaxVal = 3;
    end
    ylim([newMinVal-0.2 newMaxVal+0.2]);
    pbaspect([1.5 1 1]);
    set(gca,'fontsize',16)
    f.PaperPositionMode='auto';
    f.Position=[200 400 560 520];
    f.Color='none';
    vline(num_pos_contrasts+0.5, 'k--');


    %% create summary figures ('normalized' value of activation shift toward canonical network)

    border_data = act_mean_border(sorti, :, net_ind);
    ectopic_data = act_mean_ectopic(sorti, :, net_ind);
    canonical_data = act_canonical_mean(sorti, :, net_ind);
    noncanonical_data = act_canonical_other_mean(sorti, :, net_ind);

    % first, flip all vals where noncanon>canon
    border_data(num_pos_contrasts+1:num_contrasts, :) = border_data(num_pos_contrasts+1:num_contrasts, :) .* -1;
    ectopic_data(num_pos_contrasts+1:num_contrasts, :) = ectopic_data(num_pos_contrasts+1:num_contrasts, :) .* -1;
    canonical_data(num_pos_contrasts+1:num_contrasts, :) = canonical_data(num_pos_contrasts+1:num_contrasts, :) .* -1;
    noncanonical_data(num_pos_contrasts+1:num_contrasts, :) = noncanonical_data(num_pos_contrasts+1:num_contrasts, :) .* -1;

    % averaging values across contrasts per network
    canonical_data_mean = nanmean(canonical_data,2);
    noncanonical_data_mean = nanmean(noncanonical_data,2);
    border_data_mean = nanmean(border_data,2);
    border_norm_val = (border_data_mean - noncanonical_data_mean) ./ (canonical_data_mean - noncanonical_data_mean);
    ectopic_data_mean = nanmean(ectopic_data,2);
    ectopic_norm_val = (ectopic_data_mean - noncanonical_data_mean) ./ (canonical_data_mean - noncanonical_data_mean);    

    % create figure
    ff=figure; hold on;
    errorbar(0.9, nanmean(border_norm_val), nanstd(border_norm_val)/sqrt(num_contrasts),'bo','LineWidth',2,'markersize',5);
    errorbar(1.1, nanmean(ectopic_norm_val), nanstd(ectopic_norm_val)/sqrt(num_contrasts),'ro','LineWidth',2,'markersize',5);
    xlim([0.7 1.3])
    y_label = sprintf('norm. shift toward\ncanonical net. activ.');
    ylabel(y_label);
    title(networks{net_ind},'fontsize',26);
    pbaspect([0.4 1 1])
    set(gcf,'color','w');
    xticks([])
    yticks([0 0.5 1])
    ylim([0 1]);
    set(gca,'fontsize',18)
    ff.PaperPositionMode='auto';
    ff.Color='none';

    % save figs if desired
    if nargin==1
        saveas(f, [outDir '/HCPsubjects_' networks{net_ind} '_activations_minCanNoncanonDiff' num2str(thresh) '.fig']);
        saveas(f, [outDir '/HCPsubjects_' networks{net_ind} '_activations_minCanNoncanonDiff' num2str(thresh) '.png'], 'png');
        close(f);
        saveas(ff, [outDir '/HCPsubjects_' networks{net_ind} '_activations_SUMMARY_minCanNoncanonDiff' num2str(thresh) '.fig']);
        saveas(ff, [outDir '/HCPsubjects_' networks{net_ind} '_activations_SUMMARY_minCanNoncanonDiff' num2str(thresh) '.png'], 'png');
        close(ff);
    end
    
end