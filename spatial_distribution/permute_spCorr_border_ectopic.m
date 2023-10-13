function permute_spCorr_border_ectopic(subject_variant_info, num_perms, out_dir)

% Randomly flip a border-shift/ectopic map for a subject.
% with each permutation, half of subjects will get their maps
% randomly flipped, then overlap maps of all border-shift variants and all 
% ectopic variants are created. Spatial corr between each randomized bordershift 
% overlap map and randomized ectopic overlap map is computed.

% INPUTS
%   subject_variant_info = path to .txt file with 2 columns: subID, path to subject's border-ectopic variant map
%   num_perms = number of times to permute difference maps
%   out_dir = directory to spit out .mat with permuted difference maps


%%% read in data
[subIDs, borderEctopic] = textread(subject_variant_info, '%s%s');
num_cort_verts = 59412;


%%% load data files
all_subs_BorE = zeros(num_cort_verts, length(subIDs));
for s = 1:length(subIDs)
    BorE = ft_read_cifti_mod(borderEctopic{s});
    all_subs_BorE(:,s) = BorE.data;
end
clear borderEctopic


%%% set up permutations (flip half of subjects to switch maps)
flips = zeros(length(subIDs),1);
flips(1:(round(length(flips)/2))) = 1;


%%% loop thru permutations
spCorrs = zeros(num_perms,1);
perm_data = cell(num_perms,2);

for i=1:num_perms
    
    % initializing overlap map that will result from this permutation
    % across all subjects
    allSubsBordMap = zeros(num_cort_verts, 1);
    allSubsEcMap = zeros(num_cort_verts, 1);
    
    %decide at the SUBJECT level
    rng('shuffle'); % ensures different random seed for each permutation
    rand_flip_inds = randperm(length(flips))';
    rand_perms = flips(rand_flip_inds);
    
    for s = 1:length(subIDs)
        
        sub_borderectopic_map = all_subs_BorE(:,s);
        sub_border_map_bin = double(logical(sub_borderectopic_map==1));
        sub_ectopic_map_bin = double(logical(sub_borderectopic_map==2));

        if rand_perms(s) == 0 %do not flip
            allSubsBordMap = allSubsBordMap + sub_border_map_bin;
            allSubsEcMap = allSubsEcMap + sub_ectopic_map_bin;
        elseif rand_perms(s) == 1 %do flip
            allSubsBordMap = allSubsBordMap + sub_ectopic_map_bin;
            allSubsEcMap = allSubsEcMap + sub_border_map_bin;            
        end
    end
    
    % save the spatial corr value of overap maps for this permutation
    % also save the full permuted data
    spCorrs(i,1)= corr(allSubsBordMap, allSubsEcMap);
    perm_data{i,1}= allSubsBordMap;
    perm_data{i,2}= allSubsEcMap;
    clear allSubsBordMap allSubsEcMap
    
end

save([out_dir '/borderEctopic_spCorr_' num2str(num_perms) 'permutations.mat'],'spCorrs','perm_data')

