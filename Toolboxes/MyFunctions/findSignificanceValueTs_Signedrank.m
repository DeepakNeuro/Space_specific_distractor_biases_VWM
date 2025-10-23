function [p, p_values] = findSignificanceValueTs_Signedrank(a,b,is_paired)

% size of variables a and b should be subjects x timePts

if nargin <3
    is_paired = true;
end

trial_group_1 = a'; 
trial_group_2 = b'; 

p = nan(size(trial_group_1,1),1);
p_values = nan(size(trial_group_1,1),1);
for tt = 1:size(p,1)
    [p_values(tt), p(tt)] = signrank(trial_group_1(tt,:)', trial_group_2(tt,:)');
end
p(p==0) = nan; 

end