function [out_angles, out_N] = avg_histograms(ang_mat)

% INPUT
% ang_mat os size numSubjects x trials containing the angles within +/- 90 deg
edges_hist = -90:90; 

clear N_mat
for nSub = 1:size(ang_mat,1)
    N = histcounts(ang_mat(nSub, :), edges_hist); 
    N_mat(nSub,:) = N*(size(ang_mat, 2)/sum(N, 2)); % Normalize 
end

out_N = mean(N_mat, 1);
out_angles = [];
for e = 1:length(N)
    out_angles = [out_angles, repmat(edges_hist(e), 1, round(out_N(e)))];
end

end