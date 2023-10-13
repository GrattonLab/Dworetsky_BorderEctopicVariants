%% set up files
template_dscalar = ft_read_cifti_mod('/gpfs/research/grattonlab/HCP/Tasks/analysis/HCP_task_template.dscalar.nii');
subjects = dlmread('HCPsubjects_for_contrasts.txt'); %subjects with data for all 7 tasks and both types of variants
out_dir = '/gpfs/research/grattonlab/HCP/Tasks/analysis';
num_cortical_verts = 59412;


%% make task contrast .dscalar.nii files

task_names = {'EMOTION' 'GAMBLING' 'LANGUAGE' 'MOTOR' 'RELATIONAL' 'SOCIAL' 'WM'};
for t = 1:length(task_names)
    
    %set up task contrast info
    task = task_names{t};
    switch task
        case 'EMOTION'
            contrasts = {'FACES-SHAPES' 'FACES' 'SHAPES'};
        case 'GAMBLING'
            contrasts = {'PUNISH-REWARD' 'PUNISH' 'REWARD'};
        case 'LANGUAGE'
            contrasts = {'MATH-STORY' 'MATH' 'STORY'};
        case 'MOTOR'
            contrasts = {'LF-AVG' 'LH-AVG' 'RF-AVG' 'RH-AVG' 'T-AVG' 'AVG'};
        case 'RELATIONAL'
            contrasts = {'MATCH-REL' 'MATCH' 'REL'};
        case 'SOCIAL'
            contrasts = {'RANDOM-TOM' 'RANDOM' 'TOM'};
        case 'WM'
            contrasts = {'2BK-0BK' 'BODY-AVG' 'FACE-AVG' 'PLACE-AVG' 'TOOL-AVG' '2BK' '0BK'};
    end
    
    %compile subjects' contrast data
    all_task_subjects_data = cell(length(contrasts),1);
    all_task_subjects_mapnames = cell(length(contrasts),1);
    
    for sub = 1:length(subjects)
        subject = num2str(subjects(sub));
        sub_dir = sprintf('/gpfs/research/grattonlab/HCP/Tasks/%s/%s/MNINonLinear/Results/tfMRI_%s/tfMRI_%s_hp200_s4_level2.feat', task, subject, task, task);
        sub_cifti = sprintf('%s/%s_tfMRI_%s_level2_hp200_s4.dscalar.nii', sub_dir, subject, task);
        cifti = ft_read_cifti_mod(sub_cifti);
        
        for c = 1:length(contrasts)
            contrast = contrasts{c};
            mapname = sprintf('%s_tfMRI_%s_level2_%s_hp200_s4', subject, task, contrast);
            mapInd = find(strcmp(cifti.mapname, mapname));
            all_task_subjects_data{c}(1:num_cortical_verts, sub) = cifti.data(1:num_cortical_verts, mapInd);
            all_task_subjects_mapnames{c} = [all_task_subjects_mapnames{c} {mapname}];
        end
    end
    
    %write out contrasts as dscalar.nii files
    for c = 1:length(contrasts)
        out_data = template_dscalar;
        out_data.data = all_task_subjects_data{c};
        out_data.mapname = all_task_subjects_mapnames{c};
        out_data.hdr.dim(6) = length(subjects);
        cifti_outname = sprintf('%s/HCPsubjects_%s_%s.dscalar.nii', out_dir, task, contrasts{c});
        ft_write_cifti_mod(cifti_outname, out_data);
        clear out_data
    end
end



%% organize variant info for the above subjects

HCPsubjects_task_data_borderEctopic = zeros(num_cortical_verts, length(subjects));
HCPsubjects_task_data_reassigned = zeros(num_cortical_verts, length(subjects));
reassigned_dir = '/gpfs/research/grattonlab/HCP/Tasks/analysis/BorderEctopicData';
borderectopic_dir = '/gpfs/research/grattonlab/HCP/Tasks/analysis/BorderEctopicData';
all_task_subject_nums = cell(length(subjects), 1);

%load subjects' data
for s = 1:length(subjects)
    sub = num2str(subjects(s));
    sub_borderEctopic_file = sprintf('%s/%s_border1ectopic2.dtseries.nii', borderectopic_dir, sub);
    sub_borderEctopic = ft_read_cifti_mod(sub_borderEctopic_file);
    sub_var_networks_file = sprintf('%s/%s_reassigned.dtseries.nii', reassigned_dir, sub);
    sub_var_networks = ft_read_cifti_mod(sub_var_networks_file);
    HCPsubjects_task_data_borderEctopic(:, s) = sub_borderEctopic.data;
    HCPsubjects_task_data_reassigned(:, s) = sub_var_networks.data;
    all_task_subject_nums{s} = [all_task_subject_nums{s} sub];
    clear sub_borderEctopic sub_borderEctopic_file sub_var_networks_file sub_var_networks
end



%% write out all subjects' variant info into dscalar.nii
out_var_data = template_dscalar;
out_var_data.data = HCPsubjects_task_data_reassigned;
out_var_data.mapname = all_task_subject_nums;
out_var_data.hdr.dim(6) = length(subjects);
ft_write_cifti_mod([out_dir '/HCPsubjects_task_reassigned.dscalar.nii'], out_var_data);
out_var_data.data = HCPsubjects_task_data_borderEctopic;
ft_write_cifti_mod([out_dir '/HCPsubjects_task_borderEctopic.dscalar.nii'], out_var_data);

