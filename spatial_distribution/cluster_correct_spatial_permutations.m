function cluster_correct_spatial_permutations(perms_diffmap_mat_file, border_overlap_cifti, ectopic_overlap_cifti, num_subjects, neighbors, out_dir)

% INPUTS
%  perms_diffmap_mat_file = .mat file created by permute_diffmaps_border_ectopic.m
%  border_overlap_cifti = path to .dtseries.nii with border variant locations summed across subjects
%  ectopic_overlap_cifti = path to .dtseries.nii with ectopic variant locations summed across subjects
%  num_subjects = # of subjects used in the analysis (for thresholding)
%  neighbors = matrix containing CIFTI's neighboring vertical for cortex (assumes 7 columns)
%  out_dir = path to output CIFTIs


% setup
num_cort_verts = 59412;

% load difference maps
load(perms_diffmap_mat_file)
min_subject_threshold = ceil(0.05 * num_subjects);
diffmap_bin= logical(abs(diffmaps)>= min_subject_threshold);

% label all clusters of N>=(min_subject_threshold) subjects by diff map
all_permutations_clusters = zeros(num_cort_verts,size(diffmaps,2));

for c = 1:size(all_permutations_clusters,2)
    
    % identify discrete clusters
    clusterized = zeros(num_cort_verts,1);
    count = 1;
    for i = 1:num_cort_verts
        if diffmap_bin(i, c)

            % if the vertex is already part of a cluster, continue using
            % that ID; otherwise, use new ID
            if clusterized(i ) >0
                id = clusterized(i);
            else
                id = count;
                count = count + 1;
            end

            % label it, check its neighbors
            clusterized(i) = id;
            neighVerts = neighbors(i, 2:7);
            neighVerts(isnan(neighVerts)) = []; % Ignore NaNs

            % All nonzero neighbors are given the same ID
            for j = 1:length(neighVerts)
                % Second part of if statement prevents overwriting previously assigned variants
                if clusterized(neighVerts(j)) > 0
                    clusterized(clusterized == clusterized(neighVerts(j))) = id;
                elseif diffmap_bin(neighVerts(j), c) && clusterized(neighVerts(j)) == 0
                    clusterized(neighVerts(j)) = id;
                end
            end
        end
    end
    
    % make the IDs sequential
    ids = unique(clusterized); ids(ids == 0) = []; % Get rid of 0
    for i = 1:length(ids)
        clusterized(logical(clusterized == ids(i))) = i;
    end
      
    all_permutations_clusters(:, c) = clusterized;
    clear clusterized
end


% get all cluster sizes (in # vertices)
all_cluster_sizes = [];
for i = 1:size(all_permutations_clusters, 2)
    numClusters = length(unique(all_permutations_clusters(:, i))) - 1;
    for c = 1:numClusters
        all_cluster_sizes = [all_cluster_sizes; sum(all_permutations_clusters(:,i) == c)]; 
    end
end

% all_cluster_sizes = sort(all_cluster_sizes,'descend');
all_cluster_sizes = sort(unique(all_cluster_sizes),'descend');
top5p = round(length(all_cluster_sizes)*0.05);
cluster_threshold = all_cluster_sizes(top5p);

% read in actual difference map
b_overlap = ft_read_cifti_mod(border_overlap_cifti);
e_overlap = ft_read_cifti_mod(ectopic_overlap_cifti);
ec_minus_border = e_overlap.data - b_overlap.data;
diff_cifti = b_overlap; diff_cifti.data= ec_minus_border;
ft_write_cifti_mod([out_dir '/ectopic_minus_border_true_difference_map.dtseries.nii'], diff_cifti);

% threshold true diff map by # subjects
diff_cifti.data(abs(diff_cifti.data) < min_subject_threshold) = 0;
ft_write_cifti_mod([out_dir '/ectopic_minus_border_true_difference_map_subject_thresholded.dtseries.nii'], diff_cifti);
bin_diff_cifti= diff_cifti; bin_diff_cifti.data=logical(bin_diff_cifti.data);
ft_write_cifti_mod([out_dir '/ectopic_minus_border_true_difference_map_subject_thresholded_binarized.dtseries.nii'], bin_diff_cifti);

% re-cluster
system(['wb_command -cifti-find-clusters ' out_dir '/ectopic_minus_border_true_difference_map_subject_thresholded_binarized.dtseries.nii 0 0 0 0 COLUMN ' out_dir '/ectopic_minus_border_true_difference_map_subject_thresholded_clusters.dtseries.nii -left-surface /data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.midthickness.32k_fs_LR.surf.gii -right-surface /data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.midthickness.32k_fs_LR.surf.gii']);
true_diff_clusters = ft_read_cifti_mod([out_dir '/ectopic_minus_border_true_difference_map_subject_thresholded_clusters.dtseries.nii']);
true_diff_clustercorrected = true_diff_clusters;

% zero out clusters that do not pass size threshold
for i=1:length(unique(true_diff_clusters.data))-1
    inds = find(true_diff_clusters.data==i);
    if length(inds) < cluster_threshold
        true_diff_clustercorrected.data(inds)=0;
    end
end

true_diff_clustercorrected.data(true_diff_clustercorrected.data>0) = diff_cifti.data(true_diff_clustercorrected.data>0);

ft_write_cifti_mod([out_dir '/ectopic_minus_border_true_difference_map_subject_thresholded_cluster_corrected.dtseries.nii'], true_diff_clustercorrected);
