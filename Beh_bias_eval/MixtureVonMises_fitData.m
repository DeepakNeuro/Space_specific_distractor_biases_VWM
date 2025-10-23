%%% code gone into the papers analyiss for the von mises fit distractor
%%% alone componnent

%%
clc;
clear all;

%% Add some files to path
addpath('/media/hdd/Sanchit/Exogenous_Project/Toolboxes/MyFunctions');
addpath('/media/hdd/Sanchit/Exogenous_Project/Toolboxes/ShadedErrorBar/');
addpath('/media/hdd/Sanchit/Exogenous_Project/Data Analysis Codes/EEG Analysis/');
addpath('/media/hdd/Sanchit/Exogenous_Project/Toolboxes/bayslab_analogue_toolbox/');

%% Subjects List
sub.s2 = [108 109 110 111 113 114 115 116 117 118 119 120 121 122 123,...
    124 126 127 128 129 130 131 132 133];

Subjects =  [sub.s2];
sub.subLabel = num2str(Subjects','%d');
plt.savePath = '/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Behaviour/';

%% Initialize variables with matrices of appropriate size
for subjectNumber = Subjects
    SubIdx = find(subjectNumber==Subjects);
    fprintf('\nSubIdx = %d\n',SubIdx); 
    
    %% Load the behavioural data
    file = dir(sprintf('/media/hdd/Sanchit/Exogenous_Project/Subjects Data/Sub%d/Data/WM_S-%d_B-6*.mat',subjectNumber,subjectNumber));
    load(sprintf('%s/%s',file.folder,file.name),'trialData');
    
    %% Do Eye Tracking Rejection
    [trialData, ETrejectIdx]= doEyeTrackingRejection(subjectNumber, trialData);
    
    %% WM Response Analysis
    [angleDiff, angleDiffDist, idx_val, idx_inv, medianAbsAngleDiff, fits, actualAngle] = WMA_getAngleDiff(trialData);
    
    %% types of trials
    idx_dist_cu = (trialData.WM_cue == trialData.distractorSide); 
    idx_dist_uc = (trialData.WM_cue ~= trialData.distractorSide) & (3 ~= trialData.distractorSide);
    idx_dist_nu = (3 == trialData.distractorSide); 
    
    
    idx_interest =  ones(size(angleDiffDist)); %
    
    %% vars for selecting the trials
    idx_cu_ds = idx_val & idx_dist_cu & idx_interest;
    idx_cu_do = idx_val & idx_dist_uc & idx_interest;
    idx_cu_dn = idx_val & idx_dist_nu;
    idx_uc_ds = idx_inv & idx_dist_uc & idx_interest;
    idx_uc_do = idx_inv & idx_dist_cu & idx_interest;
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
    
    
    wm.raw.distang{1,1,SubIdx} = trialData.distractorAngle(idx_cu_ds);
    wm.raw.distang{1,2,SubIdx} = trialData.distractorAngle(idx_cu_do); 
    wm.raw.distang{1,3,SubIdx} = trialData.distractorAngle(idx_cu_dn);
    wm.raw.distang{2,1,SubIdx} = trialData.distractorAngle(idx_uc_ds); 
    wm.raw.distang{2,2,SubIdx} = trialData.distractorAngle(idx_uc_do); 
    wm.raw.distang{2,3,SubIdx} = trialData.distractorAngle(idx_uc_dn);
end


% resp angle
for sub = 1:length(Subjects)
    respang_cu_ds = [wm.raw.respang{1,1,sub}]; % Cued DS
    respang_cu_ds = respang_cu_ds(:);

    respang_cu_do = [wm.raw.respang{1,2,sub}]; % Cued DO
    respang_cu_do = respang_cu_do(:);

    respang_cu_nd = [wm.raw.respang{1,3,sub}]; % Cued no Dist
    respang_cu_nd = respang_cu_nd(:);

    % actual angle
    actang_cu_ds = [wm.raw.actang{1,1,sub}]; % Cued DS
    actang_cu_ds = actang_cu_ds(:);

    actang_cu_do = [wm.raw.actang{1,2,sub}]; % Cued DO
    actang_cu_do = actang_cu_do(:);

    actang_cu_nd = [wm.raw.actang{1,3,sub}]; % Cued no Dist
    actang_cu_nd = actang_cu_nd(:);


    % dist angle
    distang_cu_ds = [wm.raw.distang{1,1,sub}]; % Cued DS
    distang_cu_ds = distang_cu_ds(:);

    distang_cu_do = [wm.raw.distang{1,2,sub}]; % Cued DO
    distang_cu_do = distang_cu_do(:);

    distang_cu_nd = [wm.raw.distang{1,3,sub}]; % Cued no Dist
    distang_cu_nd = distang_cu_nd(:);
    
    [B_cu_ds(:,sub) LL_cu_ds(sub) ~] = mixtureFit(respang_cu_ds, actang_cu_ds);
    [B_cu_do(:,sub) LL_cu_do(sub) ~] = mixtureFit(respang_cu_do, actang_cu_do);
    [B_cu_nd(:,sub) LL_cu_nd(sub) ~] = mixtureFit(respang_cu_nd, actang_cu_nd);
    
    n_obs_cu_ds(sub) = length(~isnan(respang_cu_ds-actang_cu_ds));
    n_obs_cu_do(sub) = length(~isnan(respang_cu_do-actang_cu_do));
    n_obs_cu_nd(sub) = length(~isnan(respang_cu_nd-actang_cu_nd));
end

%% visualize the histogram

figure;
for sub_id = 1:length(Subjects)

respang_cu_ds = [wm.raw.respang{1,1,sub_id }]; % Cued DS
respang_cu_ds = respang_cu_ds;

respang_cu_do = [wm.raw.respang{1,2,sub_id }]; % Cued DO
respang_cu_do = respang_cu_do;

respang_cu_nd = [wm.raw.respang{1,3,sub_id }]; % Cued no Dist
respang_cu_nd = respang_cu_nd;

% actual angle
actang_cu_ds = [wm.raw.actang{1,1,sub_id}]; % Cued DS
actang_cu_ds = actang_cu_ds;

actang_cu_do = [wm.raw.actang{1,2,sub_id}]; % Cued DO
actang_cu_do = actang_cu_do(:);

actang_cu_nd = [wm.raw.actang{1,3,sub_id}]; % Cued no Dist
actang_cu_nd = actang_cu_nd(:);


% dist angle
distang_cu_ds = [wm.raw.distang{1,1,sub_id }]; % Cued DS
distang_cu_ds = distang_cu_ds;

distang_cu_do = [wm.raw.distang{1,2,sub_id }]; % Cued DO
distang_cu_do = distang_cu_do(:);

distang_cu_nd = [wm.raw.distang{1,3,sub_id }]; % Cued no Dist
distang_cu_nd = distang_cu_nd(:);


%%
vm_mix_fun = @(x, wU, w1, mu1, kappa1, mu2, kappa2) ...
    (wU / (2*pi)) + ...
    (w1 ./ (2*pi*besseli(0,kappa1))) .* exp(kappa1 .* cos(x - mu1)) + ...
    ((1 - wU - w1) ./ (2*pi*besseli(0,kappa2))) .* exp(kappa2 .* cos(x - mu2));


B_oi_ds =  B_cu_ds(:,sub_id);
B_oi_do =  B_cu_do(:,sub_id);
mu1 = 0;
mu2 = 2;
%% actual histogram
E = (respang_cu_ds - actang_cu_ds);
E = E(~isnan(respang_cu_ds));
E = deg2rad(E)*2;
E = mod(E + pi, 2*pi) - pi;

NE = respang_cu_ds - distang_cu_ds;
NE = NE(~isnan(respang_cu_ds));
NE = deg2rad(NE)*2;
NE = mod(NE + pi, 2*pi) - pi;
NE = NE+2;

ang = linspace(-pi, pi, 100);
kappa1 = B_oi_ds(1);kappa2 = B_oi_ds(2);wU = B_oi_ds(5); w1 = B_oi_ds(3);
y_pdf_mix3 = vm_mix_fun(ang, wU, w1, mu1, kappa1, mu2, kappa2);

subplot(4,6,sub_id);
plot(ang, y_pdf_mix3, 'LineWidth', 2); hold on;
histogram(E, linspace(-pi,pi,37), 'Normalization', 'pdf', 'DisplayStyle','stairs','EdgeColor', 'r', 'LineWidth', 0.5);
histogram(NE, linspace(-pi,pi,37), 'Normalization', 'pdf', 'DisplayStyle','stairs','EdgeColor', 'k', 'LineWidth', 0.5);
end