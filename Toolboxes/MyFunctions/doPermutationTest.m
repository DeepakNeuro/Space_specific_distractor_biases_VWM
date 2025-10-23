function [p, stats, pVal_permutationTest] = doPermutationTest(temp)

addpath('/media/hdd/Sanchit/Exogenous_Project/Toolboxes/github_repo/')

rng('default'); 
display_flag = true; % true

stats = mes(temp(:,1), temp(:,2), 'mdbysd','nBoot',10000, 'isDep', 1);

p = signrank(temp(:,1), temp(:,2),'tail','both'); %stats.t.p; 

if display_flag
fprintf('\nm1=%0.3f%c%0.3f versus m2=%0.3f%c%0.3f, P=%0.3f',...
    mean(temp(:,1)),char(177),std(temp(:,1))./sqrt(length(temp(:,1))), mean(temp(:,2)),char(177),std(temp(:,2))./sqrt(length(temp(:,2))), p); 

fprintf('\nt(%d) = %0.3f, P = %0.3f, d = %0.3f [95%% CI, %0.3f, %0.3f];',...
    stats.t.df, stats.t.tstat, stats.t.p, stats.mdbysd, stats.mdbysdCi(1), stats.mdbysdCi(2)); 
end

%% Bayes Factor
pVal_permutationTest = permutationTest(temp(:,1), temp(:,2), 1e4,'sidedness', 'both','plotresult',0);

addpath('/media/hdd/Sanchit/Exogenous_Project/Toolboxes/bayesFactor/'); 
[bf10, bf_ttest_p] = bf.ttest(temp(:,1), temp(:,2),'tail','both');
if display_flag
fprintf('\nbf10 = %0.3f; PermutationTest=%0.3f\n',...
    bf10, pVal_permutationTest); 
end
rmpath('/media/hdd/Sanchit/Exogenous_Project/Toolboxes/bayesFactor/'); 

stats.BF.bf10 = bf10; 
stats.BF.bf_ttest_p = bf_ttest_p;
end