load('permute_network_assignments_borderectopic.mat', 'permuted_ec_percentages');
num_perms = 1000;


% 1) compare true ectopic:border ratio for EACH network to the permuted
%       distribution of ratios 
% 2) assess significance via p-values (proportion of perms where true %
%       of ectopic variants is < or > (two-tailed) all perm.d percentages
%       after FDR correcting for multiple comparisons across networks

network_names = {'DMN', 'Vis', 'FP', 'DAN', 'Sal', 'CO', 'SMd', 'SMl', 'Aud', 'PMN', 'PON'};
network_IDs= [1 2 3 5 8 9 10 11 12 15 16];
networks_rgb = ...
    [0.8471, 0.1608, 0; ...
    0, 0, 0.8; ...
    0.8745, 0.8745, 0.2431; ...
    0, 0.8, 0; ...
    0, 0, 0; ...
    0.5176, 0.2039, 0.8353; ...
    0, 1, 1; ...
    1, 0.502, 0; ...
    0.7843, 0.3922, 1; ...
    0.0745, 0.3608, 0.8157; ...
    0.8, 0.8, 0.8];

%plot spread of ectopic percentages from permutations
h=figure;
hold on;
plotSpread(permuted_ec_percentages .* 100, 'distributionColors', networks_rgb, 'xnames', network_names)
set(gcf,'color','w')
ylim([0 100])
a = get(gca,'xticklabel');
set(gca,'xticklabel',a,'fontsize',16)
ylabel('Percentage of ectopic variants','fontsize',18)
set(findall(1,'type','line'),'markersize',6)
set(h, 'Position', [680 484 862 491]);

%plot true percentages on top
for i=1:length(network_IDs)
    plot(i, ec_percentages(i) .* 100, 'k.', 'MarkerSize',40)
end

%calculate network p-values based on # of random percent > or < true percent
p_values=zeros(1, length(network_IDs));
for col=1:length(network_IDs)
    true_perc = ec_percentages(col);
    permuted_percs = permuted_ec_percentages(:,col);
    num_greater = sum(permuted_percs >= true_perc);
    num_less = sum(permuted_percs <= true_perc);
    p_value= min(num_greater, num_less) / num_perms;
    p_values(col) = 2 * p_value;
end

%perform FDR correction
q = 0.05;
[~, ~, ~, adj_p] = fdr_bh(p_values, q);
sig_network_inds = find(adj_p <= q);

%plot significance (*)
for v = 1:length(sig_network_inds)
    net = sig_network_inds(v);
    plot(net,95,'*','color',[0.5 0.5 0.5], 'MarkerSize', 16);
end

pbaspect([1.8 1 1]);
