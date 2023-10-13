function [subjects_per_subgroup, subgroup_means_plot_order] = plot_variant_subgroups(assigns, corrs, meanCorr, subs_with_variants)

% Use outputs of identify_variant_subgroups.m to plot the network
% affiliation profiles of the identified subgroups. Inputs are all outputs
% from that function.
%
% OUTPUTS
%   subjects_per_subgroup = 2 columns: subjectID in first column, followed
%     by subgroupID for a given subject
%   subgroup_means_plot_order = (#subgroups) by (#networks) matrix with subgroup
%     profiles, where network order is as reordered by network_order
%     variable
%   2 figures (subgroup network profiles & subj-to-subj adjacency matrix)



% CHANGE - set network names, network order, and threshold for plotting
% (network_names should be in the same order as meanCorr (based on
% templatematch_corrCoeffs_files order)
network_names = {'DMN','Vis','FP','DAN','Sal','CO','SMd','SMl','Aud','PMN','PON'};
network_order = [9 2 7 8 6 4 3 5 1 10 11];
threshold_ind = 30;


% identify subjects in each subgroup
assigns_at_threshold = assigns(:, threshold_ind);
subgroup_IDs = unique(assigns_at_threshold);
subgroup_IDs(subgroup_IDs < 0) = [];
num_subgroups = length(subgroup_IDs);

% plot subgroups' means, sorted by network order with STE bars
% note: add more colors if num_subgroups > 8
figure;
hold on
colors = {'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'};
subgroup_means_plot_order = zeros(num_subgroups, length(network_names));
for j = 1:num_subgroups
    subgroup_ID = subgroup_IDs(j);
    group_inds = find(assigns_at_threshold == subgroup_ID);
    means2plot= meanCorr(network_order, group_inds)';
    errorbar(1:length(network_names), mean(means2plot), std(means2plot)./sqrt(size(means2plot,2)), ...
        '-','Color',colors{j},'LineWidth',2);
    subgroup_means_plot_order(j, :) = mean(meanCorr(network_order, group_inds)');
end
ylim([-0.4, 0.3]);
hline(0, 'k:')
xlim([0.5, length(network_names) + 1]);
ylabel('Mean correlation to network template', 'FontName', 'Arial', 'FontSize', 16);
set(gca, 'FontName', 'Arial', 'FontSize', 16, 'xtick', 1:length(network_names), 'xticklabel', network_names(network_order));
set(gcf,'color','w')

% plot subject-to-subject correlation matrix (grouped)
subject_corrs_order = [];
for j = 1:num_subgroups
    subgroup_ID = subgroup_IDs(j);
    subject_corrs_order = [subject_corrs_order; find(assigns_at_threshold == subgroup_ID)];
end
unassigned_sub_inds = find(assigns_at_threshold < 0);
subject_corrs_order = [subject_corrs_order; unassigned_sub_inds];

figure;
imagesc(corrs(subject_corrs_order, subject_corrs_order), [-1, 1]); 
title('Subject-to-subject variant network corrs.'); 
colormap(redblue); set(gcf,'color','w'); colorbar; axis square;


% return a variable matching subID to subgroup
subjects_per_subgroup(:, 1) = subs_with_variants;
subjects_per_subgroup(:, 2) = num2cell(assigns_at_threshold);

end

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 

end
