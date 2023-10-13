function permute_network_assignments_borderectopic(subject_variant_info, num_perms, out_dir)

% 1) permute label of each variant as border of ectopic within participants
%       (permute at participant level)
% 2) calculate proportion of pseudo-ectopic to pseudo-border variants for
%       each network
% 3) repeat (1) and (2) X times
%
% INPUTS
%   subject_variant_info = 4-column .txt file with (1) subIDs, (2) path to
%     uniqueID map, (3) path to reassigned network map, (4) path to
%     border/ectopic map
%   num_perms = # of permutations
%   out_dir = directory for .mat with permuted ectopic:border ratios/counts
%   

% load info / set-up
[subs, variantIDs, reassigned, borderectopic] = textread(subject_variant_info, '%s%s%s%s');
network_names = {'DMN' 'Vis' 'FP' 'DAN' 'Sal' 'CO' 'SMd' 'SMl' 'Aud' 'PMN' 'PON'};
network_IDs = [1 2 3 5 8 9 10 11 12 15 16];
num_cort_verts = 59412;

% pre-load data
all_subs_uids = zeros(num_cort_verts, length(subs));
all_subs_re = zeros(num_cort_verts, length(subs));
all_subs_BorE = zeros(num_cort_verts, length(subs));
for s = 1:length(subs)
    sub = subs{s};
    uids = ft_read_cifti_mod(variantIDs{s});
    re = ft_read_cifti_mod(reassigned{s});
    BorE = ft_read_cifti_mod(borderectopic{s});
    all_subs_uids(:, s) = uids.data;
    all_subs_re(:,s) = re.data;
    all_subs_BorE(:,s) = BorE.data;
end

all_perms_num_border_vars = zeros(num_perms, length(network_IDs));
all_perms_num_ectopic_vars = zeros(num_perms, length(network_IDs));

for i = 1:num_perms
    
    perm_num_border_vars = zeros(1, length(network_IDs));
    perm_num_ectopic_vars = zeros(1, length(network_IDs));
    
    for s = 1:length(subs)

        sub_uids_map = all_subs_uids(:, s);
        sub_re_map = all_subs_re(:, s);
        sub_bore_map = all_subs_BorE(:, s);
        
        num_variants = length(unique(sub_uids_map)) - 1;
        sub_net_vector = [];
        sub_borderectopic_vector = [];
        
        for v = 1:num_variants
            inds = find(sub_uids_map == v);
            sub_net_vector = [sub_net_vector; sub_re_map(inds(1))];
            sub_borderectopic_vector = [sub_borderectopic_vector; sub_bore_map(inds(1))];
        end
        
        %permute borderectopic vector
        rng('shuffle');
        rand_inds = randperm(length(sub_borderectopic_vector));
        sub_rand_borderectopic_vector = sub_borderectopic_vector(rand_inds);
        
        %keep track of randomized network/border-ectopic pairing (add across subjects)
        for n = 1:length(network_IDs)
            net_ID = network_IDs(n);
            perm_num_border_vars(1,n) = perm_num_border_vars(1,n) + sum(sub_net_vector == net_ID & sub_rand_borderectopic_vector == 1);
            perm_num_ectopic_vars(1,n) = perm_num_ectopic_vars(1,n) + sum(sub_net_vector == net_ID & sub_rand_borderectopic_vector == 2);
        end

    end
    
    all_perms_num_border_vars(i,:) = perm_num_border_vars;
    all_perms_num_ectopic_vars(i,:) = perm_num_ectopic_vars;
    
end

permuted_ec_percentages = all_perms_num_ectopic_vars ./ (all_perms_num_ectopic_vars + all_perms_num_border_vars);

save([out_dir '/permute_network_assignments_borderectopic.mat'], 'all_perms_num_border_vars', 'all_perms_num_ectopic_vars', 'permuted_ec_percentages', '-v7.3');