clc;
clear all;

%% Add some files to path
addpath('./Toolboxes/MyFunctions');
addpath('./Toolboxes/ShadedErrorBar/');
addpath('/media/hdd/Sanchit/Exogenous_Project/Data Analysis Codes/EEG Analysis/');
addpath('./Toolboxes/bayslab_analogue_toolbox/');

%% Subjects List
sub.s2 = [108 109 110 111 113 114 115 116 117 118 119 120 121 122 123,...
    124 126 127 128 129 130 131 132 133];

Subjects =  [sub.s2];

%%
for subjectNumber = Subjects
    SubIdx = find(subjectNumber==Subjects);
    fprintf('\nSubIdx = %d\n',SubIdx);
    
    %% Load the behavioural data
    file = dir(sprintf('./Subjects Data/Sub%d/Data/WM_S-%d_B-6*.mat',subjectNumber,subjectNumber));
    load(sprintf('%s/%s',file.folder,file.name),'trialData');
    
    %% Do Eye Tracking Rejection
    [trialData, ETrejectIdx]= doEyeTrackingRejection(subjectNumber, trialData);
    
    %% WM Response Analysis
    [angleDiff, idx_val, idx_inv, ~, ~, actualAngle] = WMA_getAngleDiff(trialData);
    
    %% types of trials
    idx_dist_cu = (trialData.WM_cue == trialData.distractorSide); % DS Cued side
    idx_dist_uc = (trialData.WM_cue ~= trialData.distractorSide) & (3 ~= trialData.distractorSide); % DS uncued side
    idx_dist_nu = (3 == trialData.distractorSide); % No Dist
    
    
    %% vars for selecting the trials
    idx_cu_ds = idx_val & idx_dist_cu;
    idx_cu_do = idx_val & idx_dist_uc;
    idx_cu_dn = idx_val & idx_dist_nu;
    idx_uc_ds = idx_inv & idx_dist_uc;
    idx_uc_do = idx_inv & idx_dist_cu;
    idx_uc_dn = idx_inv & idx_dist_nu;
    
    %putting in one variable
    idx_all{1,1} = idx_cu_ds;
    idx_all{1,2} = idx_cu_do;
    idx_all{1,3} = idx_cu_dn;
    idx_all{2,1} = idx_uc_ds;
    idx_all{2,2} = idx_uc_do;
    idx_all{2,3} = idx_uc_dn;
    
    % angle diff for different grouping
    wm.raw.respang{1,1,SubIdx} = trialData.ResponseWM(idx_cu_ds);
    wm.raw.respang{1,2,SubIdx} = trialData.ResponseWM(idx_cu_do);
    wm.raw.respang{1,3,SubIdx} = trialData.ResponseWM(idx_cu_dn);
    wm.raw.respang{2,1,SubIdx} = trialData.ResponseWM(idx_uc_ds);
    wm.raw.respang{2,2,SubIdx} = trialData.ResponseWM(idx_uc_do);
    wm.raw.respang{2,3,SubIdx} = trialData.ResponseWM(idx_uc_dn);
    
    
    wm.raw.actang{1,1,SubIdx} = actualAngle(idx_cu_ds);
    wm.raw.actang{1,2,SubIdx} = actualAngle(idx_cu_do);
    wm.raw.actang{1,3,SubIdx} = actualAngle(idx_cu_dn);
    wm.raw.actang{2,1,SubIdx} = actualAngle(idx_uc_ds);
    wm.raw.actang{2,2,SubIdx} = actualAngle(idx_uc_do);
    wm.raw.actang{2,3,SubIdx} = actualAngle(idx_uc_dn);
    
    wm.raw.ang{1,1,SubIdx} = angleDiff(idx_cu_ds);
    wm.raw.ang{1,2,SubIdx} = angleDiff(idx_cu_do);
    wm.raw.ang{1,3,SubIdx} = angleDiff(idx_cu_dn);
    wm.raw.ang{2,1,SubIdx} = angleDiff(idx_uc_ds);
    wm.raw.ang{2,2,SubIdx} = angleDiff(idx_uc_do);
    wm.raw.ang{2,3,SubIdx} = angleDiff(idx_uc_dn);
    
    
    wm.raw.distang{1,1,SubIdx} = trialData.distractorAngle(idx_cu_ds);
    wm.raw.distang{1,2,SubIdx} = trialData.distractorAngle(idx_cu_do);
    wm.raw.distang{1,3,SubIdx} = trialData.distractorAngle(idx_cu_dn);
    wm.raw.distang{2,1,SubIdx} = trialData.distractorAngle(idx_uc_ds);
    wm.raw.distang{2,2,SubIdx} = trialData.distractorAngle(idx_uc_do);
    wm.raw.distang{2,3,SubIdx} = trialData.distractorAngle(idx_uc_dn);
end



%% intact
%% cued
% resp angle (cued)
for subIdx = 1:24
sub_leftOut = setdiff(1:24, subIdx);
m_cu_ds = [wm.raw.ang{1,1,sub_leftOut}]; % Cued DS
m_cu_ds = m_cu_ds(:);
m_cu_ds = m_cu_ds(~isnan(m_cu_ds));

m_cu_ds = deg2rad(m_cu_ds)*2;

m_cu_do = [wm.raw.ang{1,2,sub_leftOut}]; % Cued DO
m_cu_do = m_cu_do(:);
m_cu_do = m_cu_do(~isnan(m_cu_do));
m_cu_do = deg2rad(m_cu_do)*2;

m_cu_nd = [wm.raw.ang{1,3,sub_leftOut}]; % Cued no Dist
m_cu_nd = m_cu_nd(:);
m_cu_nd = m_cu_nd(~isnan(m_cu_nd));
m_cu_nd = deg2rad(m_cu_nd)*2;

[ang, y_pdf_cu_ds, B_cu_ds(:,subIdx), wG_cu_ds] = fit_vM_simple(m_cu_ds);
[ang, y_pdf_cu_do, B_cu_do(:,subIdx), wG_cu_do] = fit_vM_simple(m_cu_do);
[ang, y_pdf_cu_nd, B_cu_nd(:,subIdx), wG_cu_nd] = fit_vM_simple(m_cu_nd);
end

%% cued req. values
delta_k_ds_nd_cu = B_cu_ds - B_cu_nd;
delta_k_ds_do_cu = B_cu_ds - B_cu_do;

MI_k_ds_nd_cu = (B_cu_ds - B_cu_nd)./(B_cu_ds + B_cu_nd);
MI_k_ds_do_cu = (B_cu_ds - B_cu_do)./(B_cu_ds + B_cu_do);

%% uncued
% resp angle (uncued)
for subIdx = 1:24
sub_leftOut = setdiff(1:24, subIdx);
m_uc_ds = [wm.raw.ang{2,1,sub_leftOut}]; % UNCued DS
m_uc_ds = m_uc_ds(:);
m_uc_ds = m_uc_ds(~isnan(m_uc_ds));

m_uc_ds = deg2rad(m_uc_ds)*2;

m_uc_do = [wm.raw.ang{2,2,sub_leftOut}]; % UNCued DO
m_uc_do = m_uc_do(:);
m_uc_do = m_uc_do(~isnan(m_uc_do));
m_uc_do = deg2rad(m_uc_do)*2;

m_uc_nd = [wm.raw.ang{2,3,sub_leftOut}]; % UNCued no Dist
m_uc_nd = m_uc_nd(:);
m_uc_nd = m_uc_nd(~isnan(m_uc_nd));
m_uc_nd = deg2rad(m_uc_nd)*2;

[ang, y_pdf_uc_ds, B_uc_ds(:,subIdx), wG_uc_ds] = fit_vM_simple(m_uc_ds);
[ang, y_pdf_uc_do, B_uc_do(:,subIdx), wG_uc_do] = fit_vM_simple(m_uc_do);
[ang, y_pdf_uc_nd, B_uc_nd(:,subIdx), wG_uc_nd] = fit_vM_simple(m_uc_nd);
end

%% uccued req. values
delta_k_ds_nd_uc = B_uc_ds - B_uc_nd;
delta_k_ds_do_uc = B_uc_ds - B_uc_do;

MI_k_ds_nd_uc = (B_uc_ds - B_uc_nd)./(B_uc_ds + B_uc_nd);
MI_k_ds_do_uc = (B_uc_ds - B_uc_do)./(B_uc_ds + B_uc_do);




%%
delta_k_ds_do_cu_mu = mean(delta_k_ds_do_cu);
delta_k_ds_do_cu_std = (sqrt(23)*std(delta_k_ds_do_cu))./sqrt(24);

delta_k_ds_do_uc_mu = mean(delta_k_ds_do_uc);
delta_k_ds_do_uc_std = (sqrt(23)*std(delta_k_ds_do_uc))./sqrt(24);



% %% Plotting code
% % Combine into arrays for plotting
% means = [delta_k_ds_do_cu_mu, delta_k_ds_do_uc_mu];
% errors = [delta_k_ds_do_cu_std, delta_k_ds_do_uc_std];
% 
% % Labels for each bar
% labels = {'DS-DO_CU', 'DS-DO_UC'};
% 
% % Create horizontal bar plot
% figure; hold on;
% b = barh(means, 'FaceColor', [0.6 0.6 0.9], 'EdgeColor', 'none');
% 
% % Add horizontal error bars
% for i = 1:numel(means)
%     errorbar(means(i), i, errors(i), 'horizontal', 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
% end
% 
% % Customize axes
% set(gca, 'YTick', 1:2, 'YTickLabel', labels, 'YDir', 'reverse'); % reverse for top-down order
% xlabel('Δk value');
% ylabel('Condition');
% title('Horizontal Bar Plot with Error Bars');
% box off;
% grid on;
% xlim padded;