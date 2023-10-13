function permute_diffmaps_border_ectopic(subject_variant_info, num_perms, out_dir)

% Randomly permute border/ectopic labels for a subject.
% with each permutation, subjects will get their border/ectopic labels
% randomly shuffled, then overlap maps of all border-shift variants and all 
% ectopic variants are created. Difference maps between each randomized 
% bordershift overlap map and randomized ectopic overlap map are created.

% INPUTS
%   subject_variant_info = path to .txt file with 3 columns: subID, path to subject's 
%       uniqueIDs variant map, path to subject's border-ectopic variant map
%   num_perms = number of times to permute difference maps
%   out_dir = directory to spit out .mat with permuted difference maps

%%% read in data
[subIDs, variantIDs, borderEctopic] = textread(subject_variant_info, '%s%s%s');
num_cort_verts = 59412;


%%% load data files
all_subs_uids = zeros(num_cort_verts, length(subs));
all_subs_BorE = zeros(num_cort_verts, length(subs));
for s = 1:length(subIDs)
    uids = ft_read_cifti_mod(variantIDs{s});
    BorE = ft_read_cifti_mod(borderEctopic{s});
    all_subs_uids(:, s) = uids.data;
    all_subs_BorE(:,s) = BorE.data;
end
clear variantIDs borderEctopic


%%% loop thru permutations
diffmaps = zeros(num_cort_verts, num_perms);

for i=1:num_perms
    
    disp(['perm ' num2str(i) '...'])
    
    %initialize permuted overlap maps
    allSubsBordMap= zeros(num_cort_verts,1);
    allSubsEcMap= zeros(num_cort_verts,1);
    
    %permute variant classifications at the subject level
    for s = 1:length(subIDs)
        
        sub_varIDs_map= all_subs_uids(:, s);
        sub_bordEc_map= all_subs_BorE(:, s);
        
        numvars = length(unique(sub_varIDs_map))-1; %variants are labeled sequentially from 1 to #vars
        rng('shuffle');
        randVarIDs = randperm(numvars);
        
        % grab the border/ectopic labels to shuffle
        bordEcLabels = zeros(numvars,1);
        for vv = 1:numvars
            varInds = find(sub_varIDs_map == vv);
            bordEcLabels(vv,1) = sub_bordEc_map(varInds(1));
        end
        
        % shuffle them by random ID
        randBordEcLabels = zeros(numvars, 1);
        randBordMap = zeros(num_cort_verts, 1);
        randEcMap = zeros(num_cort_verts, 1);
        for v = 1:numvars
            randBordEcLabels(v,1) = bordEcLabels(randVarIDs(v));
            
            % create new permuted subject-level map
            if randBordEcLabels(v,1) == 1
                randBordMap(sub_varIDs_map == v) = 1;
            elseif randBordEcLabels(v,1) == 2
                randEcMap(sub_varIDs_map == v) = 1;
            end
        end
        
        % add subject's map to overall permutation overlap map
        allSubsBordMap= allSubsBordMap + randBordMap;
        allSubsEcMap = allSubsEcMap + randEcMap;
    end
    
    % save this permutation's difference map
    diffmaps(:,i) = allSubsEcMap - allSubsBordMap;
    
end

save([out_dir '/borderEctopic_diffmaps_' num2str(num_perms) 'permutations.mat'], 'diffmaps')