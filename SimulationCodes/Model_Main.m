clc; clear all; 

%% Add to path
addpath('/media/hdd/Sanchit/Exogenous_Project/Toolboxes/MyFunctions');
addpath('/media/hdd/Sanchit/Exogenous_Project/Data Analysis Codes/');
addpath('/media/hdd/Sanchit/Exogenous_Project/Data Analysis Codes/EEG Analysis/');
plt.savePath = '/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Simulations_DistractorProject/';

rng('default')
nIter = 300; %1e4;  %2*1e4;%1*1e3;

clear sim; 
sim.ori_mat = nan(nIter,3);
sim.dist_params_mat = nan(nIter,2);
sim.ang_cu_mat = nan(nIter,1525,2);
sim.ang_uc_mat = nan(nIter,1525,2);


rng('default'); 

dist_time_distribution = exprnd(5, 1500,1)*(1000/15);
dist_time_distribution = 1700 + dist_time_distribution(dist_time_distribution<=750);

for iter = 1:nIter
    disp(iter);
    rng(iter);
    ori = randsample([-90:90], 3, true);
    dist_time = randsample(dist_time_distribution, 1);
    dist_side = randsample([1, 2, 3], 1);
    % simulation main, change this funciton to simulate different parameter
    % values for encoding, maintenence, dist. encoding, dist. gating
    [ang_cu, ang_uc, metrics] = simulation_function(ori, [dist_time, dist_side], iter);
    
    sim.ori_mat(iter,:) = ori;
    sim.dist_params_mat(iter,:) = [dist_time, dist_side];
    sim.ang_cu_mat(iter,:,:) = ang_cu;
    sim.ang_uc_mat(iter,:,:) = ang_uc;
end

time = 0:2:3050-2; 

save(sprintf('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Simulations_DistractorProject/Simulations_%s.mat',dateTime), '-v7.3');