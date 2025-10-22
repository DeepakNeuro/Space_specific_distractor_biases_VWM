function [p, p_values, t_sums] = findSignificanceValueTs_ClusterPermutest(a,b,is_paired)

% size of variables a and b should be subjects x timePts

rng('default');

if nargin <3
    is_paired = true;
end

trial_group_1 = a'; 
trial_group_2 = b'; 

[clusters, p_values, t_sums, permutation_distribution ] = permutest(trial_group_1 , trial_group_2, is_paired, ...
    0.05, 10^4, false);
p = nan(size(trial_group_1,1),1);
for kk = 1:size(clusters,2)
    if p_values(kk) <= p_values(1)
        p(clusters{kk}) = 1;
    end
end

end