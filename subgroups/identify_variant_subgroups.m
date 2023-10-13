function [meanCorr, corrs, assigns, subs_with_variants] = identify_variant_subgroups(sublist, templatematch_corrCoeffs, out_dir)

% Use Infomap to subgroup individuals based on network affiliations of
% variant locations.
% Intended to be run after SeitzmanGratton-2019-PNAS script templateMatchingVariants
% script, separately for border and ectopic variants; see
% https://github.com/GrattonLab/SeitzmanGratton-2019-PNAS/blob/master/templateMatchingVariants.m
% 
% INPUTS:
%   sublist = .txt list of subject IDs (must match #/order of maps in each CIFTI
%     in templatematch_corrCoeffs files below)
% 
%   templatematch_corrCoeffs = .txt file with each line pointing to a
%   network corrCoeff file (each with #subjects maps); e.g.:
%       /path/to/Templatematch_corrCoeffs_borderOnly_DMN.dtseries.nii
%       /path/to/Templatematch_corrCoeffs_borderOnly_Vis.dtseries.nii
%       /path/to/Templatematch_corrCoeffs_borderOnly_FP.dtseries.nii
%       /path/to/Templatematch_corrCoeffs_borderOnly_Sal.dtseries.nii
%       /path/to/Templatematch_corrCoeffs_borderOnly_CO.dtseries.nii
%       ...
%   out_dir = directory for infomap outputs
%
% OUTPUTS
%   meanCorr = contains network similarity vector for each subject
%   corrs = contains subject-to-subject adjacency matrix of meanCorr
%     (infomap input)
%   assigns = infomap-derived subgroup assignments across thresholds
%     specified in script
%   goodsubs = sublist, excluding subjects with no variants, if applicable
%     (e.g., if a given subject has no ectopic variants). meanCorr and corrs
%     will also reflect these subject exclusions
%


% read file inputs
templatematch_corrCoeffs_files = textread(templatematch_corrCoeffs, '%s');
subs = textread(sublist, '%s');
subs_with_variants = subs;


% compute subjects' mean correlation coeff. to each template network 
% (variant-size-weighted; across variant vertices)
meanCorr = zeros(length(templatematch_corrCoeffs_files), length(subs));
for i = 1:length(templatematch_corrCoeffs_files)
    fileName = templatematch_corrCoeffs_files{i};
    dataTemp = ft_read_cifti_mod(fileName);
    dataTemp = dataTemp.data;
    dataTemp(dataTemp==0) = NaN;
    meanCorr(i,:) = nanmean(dataTemp);
end


% remove subs with no variants (same subj. indices for all corrCoeff files)
bad_subs = find(isnan(meanCorr(1,:)));
subs_with_variants(bad_subs)  =[];
meanCorr(:, bad_subs) = [];
corrs = corr(meanCorr);


% change Infomap parameters if desired
threshold_array = 0.05:0.01:0.5;
min_subgroup_size = 20;  %minimum number of subjects required for a subgroup


% run Infomap to subgroup
cd(out_dir);
Run_Infomap_2015(corrs, ones(size(corrs)), 0, threshold_array, 0, out_dir)
modify_clrfile('simplify','rawassn.txt',min_subgroup_size); %makes rawassn_minsizeX.txt
assigns = dlmread(['rawassn_minsize' num2str(min_subgroup_size) '.txt']);


% plot results across thresholds
figure; imagesc(assigns); 
title('Infomap assignments across thresholds');
set(gcf, 'color', 'w'); 
set(gca, 'XTick', 5:5:length(threshold_array))
set(gca, 'XTickLabel', threshold_array(5:5:end));
ylabel('subjects');
xlabel('thresholds');
saveas(gcf, [out_dir '/infomap_assignments'], 'fig')
saveas(gcf, [out_dir '/infomap_assignments'], 'epsc')
close(gcf);

end


