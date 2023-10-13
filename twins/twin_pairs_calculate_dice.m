%% load subj. info and variant info
MZsubsInTwins = dlmread('MZsubsInTwins.txt'); %2 columns of paired subIDs
DZsubsInTwins = dlmread('DZsubsInTwins.txt'); %2 columns of paired subIDs
nontwin_pairs = dlmread('nontwin_pairs.txt'); %2 columns of paired subIDs
remaining_subs_unrelated_paired = dlmread('remaining_subs_unrelated_paired.txt'); %2 columns of paired subIDs
borderEctopicLoc = '/data/cn6/allyd/BorderEctopic/twins_fullcohort/borderEctopic';
parcFree = false;

%% calculate Dice for each pairing
if ~parcFree
    add_str = '';
elseif parcFree
    add_str = '_parcFree';
end

disp('MZ')
MZ_border_dice=zeros(size(MZsubsInTwins,1),1);
MZ_ectopic_dice=zeros(size(MZsubsInTwins,1),1);
for i=1:size(MZsubsInTwins,1)
    sub1=MZsubsInTwins(i,1); sub2=MZsubsInTwins(i,2);
    sub1_be= ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub1) '_border1ectopic2' add_str '.dtseries.nii']); 
    sub2_be= ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub2) '_border1ectopic2' add_str '.dtseries.nii']); 
    MZ_border_dice(i,1)= getDiceCoeff(sub1_be.data==1, sub2_be.data==1);
    MZ_ectopic_dice(i,1)= getDiceCoeff(sub1_be.data==2, sub2_be.data==2);
end

disp('DZ')
DZ_border_dice=zeros(size(DZsubsInTwins,1),1);
DZ_ectopic_dice=zeros(size(DZsubsInTwins,1),1);
for i=1:size(DZsubsInTwins,1)
    sub1=DZsubsInTwins(i,1); sub2=DZsubsInTwins(i,2);
    sub1_be= ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub1) '_border1ectopic2' add_str '.dtseries.nii']); 
    sub2_be= ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub2) '_border1ectopic2' add_str '.dtseries.nii']); 
    DZ_border_dice(i,1)= getDiceCoeff(sub1_be.data==1, sub2_be.data==1);
    DZ_ectopic_dice(i,1)= getDiceCoeff(sub1_be.data==2, sub2_be.data==2);
end

disp('SIBS')
SIBS_border_dice=zeros(size(nontwin_pairs,1),1);
SIBS_ectopic_dice=zeros(size(nontwin_pairs,1),1);
for i=1:size(nontwin_pairs,1)
    sub1=nontwin_pairs(i,1); sub2=nontwin_pairs(i,2);
    sub1_be= ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub1) '_border1ectopic2' add_str '.dtseries.nii']); 
    sub2_be= ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub2) '_border1ectopic2' add_str '.dtseries.nii']); 
    SIBS_border_dice(i,1)= getDiceCoeff(sub1_be.data==1, sub2_be.data==1);
    SIBS_ectopic_dice(i,1)= getDiceCoeff(sub1_be.data==2, sub2_be.data==2);
end

disp('UNR')
UNR_border_dice=zeros(size(remaining_subs_unrelated_paired,1),1);
UNR_ectopic_dice=zeros(size(remaining_subs_unrelated_paired,1),1);
for i=1:size(remaining_subs_unrelated_paired,1)
    sub1=remaining_subs_unrelated_paired(i,1); sub2=remaining_subs_unrelated_paired(i,2);
    sub1_be= ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub1) '_border1ectopic2' add_str '.dtseries.nii']); 
    sub2_be= ft_read_cifti_mod([borderEctopicLoc '/' num2str(sub2) '_border1ectopic2' add_str '.dtseries.nii']); 
    UNR_border_dice(i,1)= getDiceCoeff(sub1_be.data==1, sub2_be.data==1);
    UNR_ectopic_dice(i,1)= getDiceCoeff(sub1_be.data==2, sub2_be.data==2);
end