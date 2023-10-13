%% description of variables

% MZorDZ: (#pairs)x1 matrix indicating whether the pair is an MZ or DZ pair (MZ is 0; DZ is 1)
% border_varmaps: (#subjects)x2 cell array; each cell is a 59412x1 binarized
%    map of border variants for the corresponding subject in subsInTwins
% ectopic_varmaps: same as above, for ectopic variants


outdir = '/data/cn6/allyd/BorderEctopic/removeBlueBrains/heritability/permutations';
numperms = 1000;


%load subject pair lists (2 columns with each row consisting of 2 paired subIDs)
MZsubsInTwins = dlmread('MZsubsInTwins.txt');
DZsubsInTwins = dlmread('DZsubsInTwins.txt');
borderEctopicLoc = '/data/cn6/allyd/BorderEctopic/twins_fullcohort/borderEctopic';
parcFree = false;

%% set up for permutations

%make MZ/DZ labels (0=MZ, 1=DZ)
subsInTwins = [MZsubsInTwins; DZsubsInTwins];
MZorDZ = zeros(size(subsInTwins,1),1);
MZorDZ(size(MZsubsInTwins,1)+1:end) = 1; 

%load corresponding variant maps and calculate true dice vals between pairs
border_varmaps = cell(size(subsInTwins,1),2);
ectopic_varmaps = cell(size(subsInTwins,1),2);
true_border_dice = [];
true_ectopic_dice = [];

for i=1:size(subsInTwins,1)
    sub1 = subsInTwins(i,1);
    sub2 = subsInTwins(i,2);
    if ~parcFree
        sub1be = ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub1) '_border1ectopic2.dtseries.nii']);
        sub2be = ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub2) '_border1ectopic2.dtseries.nii']);
    elseif parcFree
        sub1be = ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub1) '_border1ectopic2_parcFree.dtseries.nii']);
        sub2be = ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub2) '_border1ectopic2_parcFree.dtseries.nii']);
    end
    border_varmaps{i,1} = double(sub1be.data==1);
    ectopic_varmaps{i,1} = double(sub1be.data==2);
    border_varmaps{i,2} = double(sub2be.data==1);
    ectopic_varmaps{i,2} = double(sub2be.data==2);
    true_border_dice = [true_border_dice; getDiceCoeff(border_varmaps{i,1}, border_varmaps{i,2})];
    true_ectopic_dice = [true_ectopic_dice; getDiceCoeff(ectopic_varmaps{i,1}, ectopic_varmaps{i,2})];
end


%% permute for border variants and ectopic variants separately

vartypes = {'border', 'ectopic'}; 

for v = 1:2
    if strcmp(vartypes{v}, 'border')
        varmaps = border_varmaps;
        perm_dice_vals = true_border_dice;
    elseif strcmp(vartypes{v}, 'ectopic')
        varmaps = ectopic_varmaps;
        perm_dice_vals = true_ectopic_dice;
    end

    %begin permutation of MZ/DZ labels and recalculate Falconer's 
    perm_falconers = zeros(numperms,1);
    perm_all_dice = zeros(size(subsInTwins,1), numperms);
    perm_all_labels = zeros(size(subsInTwins,1), numperms);
    for i = 1:numperms
        rng('shuffle');
        iperm = randperm(size(subsInTwins,1))';

        perm_labels = MZorDZ(iperm);

        %permute pairings
        MZ_perm_dice_vals = perm_dice_vals(perm_labels==0);
        DZ_perm_dice_vals = perm_dice_vals(perm_labels==1);

        perm_falconers(i) = 2 * (mean(MZ_perm_dice_vals) - mean(DZ_perm_dice_vals));
        perm_all_labels(:,i) = perm_labels;
    end
    save([outdir '/' vartypes{v} '_' num2str(numperms) 'perms.mat'], 'perm_falconers', 'perm_all_labels');
end


%% calculate true Falconers values
falconers_border = 2 * (mean(true_border_dice(MZorDZ==0)) - mean(true_border_dice(MZorDZ==1)));
falconers_ectopic = 2 * (mean(true_ectopic_dice(MZorDZ==0)) - mean(true_ectopic_dice(MZorDZ==1)));


%% plot
load([outdir '/border_' num2str(numperms) 'perms.mat'], 'perm_falconers');
perm_falconers_border = perm_falconers; clear perm_falconers
load([outdir '/ectopic_' num2str(numperms) 'perms.mat'], 'perm_falconers');
perm_falconers_ectopic = perm_falconers; clear perm_falconers

figure; hold on;
plotSpread(perm_falconers_border,'distributionColors',[0.5 0.5 0.5],'markersize',15)
plotSpread(falconers_border,'distributionColors','r','markersize',40)
set(gcf,'color','w')
set(gca,'fontsize',18,'fontweight','bold')
set(gca,'xTick',[])
xlim([0.25 1.75])
ylim([-0.15002 0.15002])
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0 0 8 15])
title({'Permuted Falconer''s:','Border variant overlap'},'FontSize',20,'fontweight','normal')
saveas(gcf, [outdir '/border_' num2str(numperms) 'perms'], 'fig')
saveas(gcf, [outdir '/border_' num2str(numperms) 'perms'], 'png')
close(gcf)

figure; hold on;
plotSpread(perm_falconers_ectopic,'distributionColors',[0.5 0.5 0.5],'markersize',15)
plotSpread(falconers_ectopic,'distributionColors','r','markersize',40)
set(gcf,'color','w')
set(gca,'fontsize',18,'fontweight','bold')
set(gca,'xTick',[])
xlim([0.25 1.75])
ylim([-0.15002 0.15002])
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[0 0 8 15])
title({'Permuted Falconer''s:','Ectopic variant overlap'},'FontSize',20,'fontweight','normal')
saveas(gcf, [outdir '/ectopic_' num2str(numperms) 'perms'], 'fig')
saveas(gcf, [outdir '/ectopic_' num2str(numperms) 'perms'], 'png')
close(gcf)


