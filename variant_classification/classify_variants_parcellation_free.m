function classify_variants_parcellation_free(subject_data_info, group_dconn_loc, dists, outdir)

%
% This function classifies a variant as border or ectopic based on its
% seedmap's peak R-to-group-seedmap values as the distance increases between the 
% variant and vertices in the group-average. The codebase for producing
% necessary inputs (variants labeled by ID and by network) can be found at:
% https://github.com/GrattonLab/SeitzmanGratton-2019-PNAS
% This script requires supporting scripts for reading and writing CIFTI
% files, where the modified versions contained in this script can be found
% at the above link.
%
%
% INPUTS
% subject_data_info = path to a .txt file with 3 columns, each row containing: subject, uniqueIDs variant map location, dconn map location; e.g.:
%     MSC01    /path/to/variantIDs/MSC01_uniqueIDs_afterReassign.dtseries.nii    /path/to/dconn/MSC01_corrs.dconn.nii
%     MSC02    /path/to/variantIDs/MSC02_uniqueIDs_afterReassign.dtseries.nii    /path/to/dconn/MSC02_corrs.dconn.nii
%       ...
% group_dconn_loc = path to the dconn file (vertex-to-vertex connectivity matrix CIFTI) of the
%     group average (e.g., '/path/to/groupavg/120_avg_corr_LR.dconn.nii')
% dists = array of distances at which to calculate variant's peak R to group-avg. locations (e.g., [10 150] - must include 10 and 150)
% outdir = directory for output .mat file
%
% OUTPUTS
% {sub}_variants_peak_corr_to_group.mat, containing:
%   variant_peaks_cumulative max = max variant-R-to-group at each distance (rows=variants; columns=distances)
%   dists = original array of distances input
%   border_variants_parc_free = vector with variant uniqueIDs corresponding to border variants
%   ectopic_variants_parc_free = vector with variant uniqueIDs corresponding to ectopic variants
% a CIFTI file for each subject where border variant locations are labeled
%   with a value of 1 and ectopic variant locations are labeled with a
%   value of 2
%

%% setup

num_cort_verts = 59412; %cortical vertices for resizing CIFTIs

% check dists input
if ~ismember(10, dists) || ~ismember(150, dists)
    error('10 and 150 must be included in ''dists'' variable to perform parcellation-free classification')
end
% load/resize distance matrix
disp('loading distance matrix...')
dmat_loc = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_geo_distances_xhemisphere_large.mat';
load(dmat_loc); %loads [distances] variable (vertex-to-vertex)
distances= single(distances(1:num_cort_verts, 1:num_cort_verts));


%% load subject information and begin loop

[subs, variantIDs_loc, dconn_loc] = textread(subject_data_info, '%s%s%s');

% load group average correlation matrix
disp('loading group avg. dconn...')
group_dconn = ft_read_cifti_mod(group_dconn_loc);
group_dconn = group_dconn.data(1:num_cort_verts, 1:num_cort_verts);
disp('group avg. dconn loaded')

for s = 1:length(subs)
    
    sub = subs{s};
    disp([sub ': classifying variants (parcellation-free)'])
    
    % load dconn and variant IDs map
    disp('  loading subject dconn...')
    dconn = ft_read_cifti_mod(dconn_loc{s});
    dconn = dconn.data(1:num_cort_verts, 1:num_cort_verts);
    disp('  subject dconn loaded')
    var_ids = ft_read_cifti_mod(variantIDs_loc{s});
    b1e2 = var_ids; b1e2.data = zeros(size(var_ids.data,1), 1);
    
    
    variant_peaks = zeros(length(unique(var_ids.data)) - 1, length(dists)); 
    
    % loop through variants
    for var = 1:length(unique(var_ids.data)) - 1
        inds = find(var_ids.data == var);

        layers_inds = cell(length(dists), 1);
        measured_inds = inds;
 

        % calculate variant's average seedmap
        var_avg_seedmap = mean(dconn(inds,:), 1);

        % keep track of vertices within X distances from variant (exclusive
        % of each other)
        for d = 1:length(dists)
            dist = dists(d);
            
            close_verts = [];
            for ind = 1:length(inds)
                close_verts = [close_verts; find(distances(inds(ind), :) <= dists(d))'];
            end

            outer_layer= setdiff(close_verts, measured_inds);
            layers_inds{d} = outer_layer;

            if dist == 0 %variant location only - no surrounding vertices
                layers_inds{d} = inds;
            end
            measured_inds=[measured_inds; outer_layer];
        end

        peaks = [];
        for d = 1:length(dists)
            corrs = zeros(length(layers_inds{d}), 1);
            for oind = 1:length(layers_inds{d})
                nInd = layers_inds{d}(oind);
                corrs(oind, 1) = corr(var_avg_seedmap', group_dconn(nInd, :)');
            end
            peaks = [peaks; max(corrs)]; 
        end
        variant_peaks(var, :) = peaks;
    end
    
    % calculate cumulative max across distances
    variant_peaks_cumulative_max = cummax(variant_peaks, 2);
    
    
    % loop through variants again and classify as ectopic if R at 10mm is
    % >=90% of peak (as calculated by max R at 150mm)
    criteria_ind = find(dists == 10);
    max_dist_ind = find(dists == 150);

    border_variants_parc_free = [];
    ectopic_variants_parc_free = [];
    for vv = 1:size(variant_peaks_cumulative_max, 1)
        min_val = 0.9 * variant_peaks_cumulative_max(vv, max_dist_ind);
        if variant_peaks_cumulative_max(vv, criteria_ind) >= min_val
            border_variants_parc_free = [border_variants_parc_free; vv];
            b1e2.data(var_ids.data == vv) = 1;
        elseif variant_peaks_cumulative_max(vv, criteria_ind) < min_val
            ectopic_variants_parc_free = [ectopic_variants_parc_free; vv];
            b1e2.data(var_ids.data == vv) = 2;
        end
    end
    
    % save peak values at distances and new border/ectopic variant indices
    save([outdir '/' sub '_variants_peak_corr_to_group.mat'], 'dists', 'variant_peaks_cumulative_max', 'border_variants_parc_free', 'ectopic_variants_parc_free', '-v7.3');
    ft_write_cifti_mod([outdir '/' sub  '_border1ectopic2_parcFree.dtseries.nii'], b1e2)

   
end

end