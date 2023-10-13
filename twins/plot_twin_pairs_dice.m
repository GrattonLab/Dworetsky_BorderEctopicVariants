function plot_twin_pairs_dice(MZ_border_dice, MZ_ectopic_dice, DZ_border_dice, DZ_ectopic_dice, SIBS_border_dice, SIBS_ectopic_dice, UNR_border_dice, UNR_ectopic_dice)

%% organize data to plot

border_pairs{1}=MZ_border_dice; border_pairs{2}=DZ_border_dice; border_pairs{3}=SIBS_border_dice; border_pairs{4}=UNR_border_dice;
ectopic_pairs{1}=MZ_ectopic_dice; ectopic_pairs{2}=DZ_ectopic_dice; ectopic_pairs{3}=SIBS_ectopic_dice; ectopic_pairs{4}=UNR_ectopic_dice;

group=[ones(1, length(MZ_border_dice)), ones(1, length(MZ_ectopic_dice)).*2, ones(1, length(DZ_border_dice)).*3, ...
    ones(1, length(DZ_ectopic_dice)).*4, ones(1, length(SIBS_border_dice)).*5, ones(1, length(SIBS_ectopic_dice)).*6, ...
    ones(1, length(UNR_border_dice)).*7 ones(1, length(UNR_ectopic_dice)).*8];
data= [border_pairs{1}' ectopic_pairs{1}' border_pairs{2}' ectopic_pairs{2}' border_pairs{3}' ectopic_pairs{3}' border_pairs{4}' ectopic_pairs{4}'];

MZ_border_sem = std(MZ_border_dice) ./ sqrt(length(MZ_border_dice));
MZ_ectopic_sem = std(MZ_ectopic_dice) ./ sqrt(length(MZ_ectopic_dice));
DZ_border_sem = std(DZ_border_dice) ./ sqrt(length(DZ_border_dice));
DZ_ectopic_sem = std(DZ_ectopic_dice) ./ sqrt(length(DZ_ectopic_dice));
SIBS_border_sem = std(SIBS_border_dice) ./ sqrt(length(SIBS_border_dice));
SIBS_ectopic_sem = std(SIBS_ectopic_dice) ./ sqrt(length(SIBS_ectopic_dice));
UNR_border_sem = std(UNR_border_dice) ./ sqrt(length(UNR_border_dice));
UNR_ectopic_sem = std(UNR_ectopic_dice) ./ sqrt(length(UNR_ectopic_dice));

MZ_DZ_SIBS_UNR_border_means = [mean(MZ_border_dice) mean(DZ_border_dice) mean(SIBS_border_dice) mean(UNR_border_dice)];
MZ_DZ_SIBS_UNR_border_sems = [MZ_border_sem DZ_border_sem SIBS_border_sem UNR_border_sem];
MZ_DZ_SIBS_UNR_ectopic_means = [mean(MZ_ectopic_dice) mean(DZ_ectopic_dice) mean(SIBS_ectopic_dice) mean(UNR_ectopic_dice)];
MZ_DZ_SIBS_UNR_ectopic_sems = [MZ_ectopic_sem DZ_ectopic_sem SIBS_ectopic_sem UNR_ectopic_sem];
positions=[1 1.25 2 2.25 3 3.25 4 4.25];

figure;errorbar([1:4]-0.05,MZ_DZ_SIBS_UNR_border_means,MZ_DZ_SIBS_UNR_border_sems,'bo','LineWidth',2); hold on;
errorbar([1:4]+0.05,MZ_DZ_SIBS_UNR_ectopic_means,MZ_DZ_SIBS_UNR_ectopic_sems,'ro','LineWidth',2); hold on;
title({'Subject pairs: Dice of variant locations','Network-dependent classification'},'fontsize',16)
set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) mean(positions(5:6)) mean(positions(7:8))],'fontsize',14)
set(gca,'xticklabel',{'MZ','DZ','SIBS','Unrelated'})
ylabel('subject-to-subject dice','fontsize',16)
ylim([0 0.2]); set(gcf,'color','w')
xlim([0.75 4.5])

end