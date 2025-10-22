%{
Code for reproducing Fig. 2 of the paper. 
%}
Figures_main;
%% Load EEG data
% Load the memorandum decoding output
eeg.dec.stim = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/EEG_posterior_all/AngSpaces_longPreStim_Stim_3_6s_epoch_badelecs_interpl_longest_HemiWise_39elecs_EEGAnalysis_Decoding_June_17_2024_15_25_34');

% x axis of tuning curve
step_size = 180/16;
angspace_ang2 = (90:-step_size:-90) - step_size/2;
angspace_ang2(end) = [];

%% Fig. 2A. Decoding tuning map as a function of time 
time = eeg.dec.stim.eeg.decodingTrC.StimulusOnset.time(1,:);
dis_mem_l = mean(nanmean(eeg.dec.stim.eeg.decodingTrC.StimulusOnset.l.er.dists_mem,5),3); 
dis_mem_r = mean(nanmean(eeg.dec.stim.eeg.decodingTrC.StimulusOnset.r.el.dists_mem,5),3);

dis_mem = squeeze(mean((dis_mem_l + dis_mem_r)/2,1));

% figure decoing tuning profile map
color_map_dec = flip(brewermap(250,'RdBu'));
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact'); 
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [2 3];

tp = nexttile(tl1, 1);
imagesc(time, angspace_ang2, dis_mem); box off; 
colormap(color_map_dec); cb = colorbar; 
caxis([-6 6]*1e-3);
xlabel('Time from stimulus onset (s)')
ylabel('Angle Difference (°)')
hold on;

xlim([-0.1, 0.57])
ylim([-1, 1]*90)
ax = gca; 
ax.Box = 'off';
ax.XTick = [0:0.25:0.5];
ax.YTick = [-90:45:90];

%% Fig. 2B. Memorandum decoding strength 
delayDuration = beh.trialData.delayDuration+1.2;
time = eeg.dec.stim.eeg.decodingTrC.StimulusOnset.time(1,:);
idxL = eeg.dec.stim.eeg.decodingTrC.beh.attendL(1,:);
idxR = eeg.dec.stim.eeg.decodingTrC.beh.attendR(1,:);
idxL = idxL;
idxR = idxR;
%%
cos_mem_l = mean(nanmean(eeg.dec.stim.eeg.decodingTrC.StimulusOnset.l.er.cos_mem, 4),2);
cos_mem_r = mean(nanmean(eeg.dec.stim.eeg.decodingTrC.StimulusOnset.r.el.cos_mem, 4),2);
cos_mem_avg = squeeze((cos_mem_l + cos_mem_r)/2);
%%
cos_mem_l_al = squeeze(mean(nanmean(eeg.dec.stim.eeg.decodingTrC.StimulusOnset.l.er.cos_mem(:,:,:,idxL),4),2));
cos_mem_r_al = squeeze(mean(nanmean(eeg.dec.stim.eeg.decodingTrC.StimulusOnset.r.el.cos_mem(:,:,:,idxL),4),2));
cos_mem_l_ar = squeeze(mean(nanmean(eeg.dec.stim.eeg.decodingTrC.StimulusOnset.l.er.cos_mem(:,:,:,idxR),4),2));
cos_mem_r_ar = squeeze(mean(nanmean(eeg.dec.stim.eeg.decodingTrC.StimulusOnset.r.el.cos_mem(:,:,:,idxR),4),2));

cos_mem_avg_cu = (cos_mem_l_al + cos_mem_r_ar)/2;
cos_mem_avg_uc = (cos_mem_l_ar + cos_mem_r_al)/2;

cos_mem = smoothdata(mean(cos_mem_avg,1),'movmean', 10);
cos_mem_se = std(cos_mem_avg,[],1)./sqrt(nEEG);

% cluster based permutaion test
[pVal_line, pVal] = findSignificanceValueTs_ClusterPermutest(cos_mem_avg, cos_mem_avg*0, true);

cos_mem_cu = smoothdata(mean(cos_mem_avg_cu,1),'movmean', 10);
cos_mem_uc = smoothdata(mean(cos_mem_avg_uc,1),'movmean', 10);
cos_mem_cu_se = std(cos_mem_avg_cu,[],1)./sqrt(nEEG);
cos_mem_uc_se = std(cos_mem_avg_uc,[],1)./sqrt(nEEG);

% [pVal_line_cu, pVal_cu] = findSignificanceValueTs_ClusterPermutest(cos_mem_avg_cu, cos_mem_avg_cu*0, true);
% [pVal_line_uc, pVal_uc] = findSignificanceValueTs_ClusterPermutest(cos_mem_avg_uc, cos_mem_avg_uc*0, true);

% figures
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact'); 
tl1.Layout.Tile = 13;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1,1);

hold on; 
ylim([-4, 15]*1e-4)
xlim([-0.1, 0.57])
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = boundedline(time, cos_mem, cos_mem_se, 'cmap', color_map(6,:));
plot(time, pVal_line*(-2*1e-4), 'Color', color_map(6,:), 'LineWidth', 1);
ax = gca; 
ax.Box = 'off';
ax.XTick = [0:0.25:0.5];
ax.YTick = [0 10]*1e-4;
xlabel('Time (s)'); 
title('Stimulus Decoding')

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact'); 
tl1.Layout.Tile = 15;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1,1);
hold on; 
ylim([-2, 4.5]*1e-4)
xlim([0.1, 1.2])
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;


p = boundedline(time, [cos_mem_cu; cos_mem_uc]', [cos_mem_cu_se; cos_mem_uc_se]', 'alpha','cmap', color_map_att); hold on;
plot(time, pVal_line_cu*(-1*1e-4), 'Color', color_map(6,:), 'LineWidth', 1);
plot(time, pVal_line_uc*(-1.2*1e-4), 'Color', color_map(6,:), 'LineWidth', 1);
ax = gca; 
ax.Box = 'off';
xlabel('Time (s) from stim onset'); 
title('Stimulus Decoding')
ylabel('Decoding Accuracy (norm. µV)', 'FontSize', 8); 



%% violin check for cued Vs uncued  
ts_win = time>=0.7 & time<=0.8;

cos_mem_cu_toi = nanmean(cos_mem_avg_cu(:, ts_win),2);
cos_mem_uc_toi = nanmean(cos_mem_avg_uc(:, ts_win),2);

m = [cos_mem_uc_toi,zeros(size(cos_mem_uc_toi))];
doPermutationTest(m,'larger','right');

% figure
tl1 = tiledlayout(tl, 1,1, 'TileSpacing','compact', 'Padding', 'compact'); 
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [3 2];

tp = nexttile(tl1,1);

xlim([0.5 2.5])
ylim([-3 3]*1e-3);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
hold on;
[vp, bp] = violin_plot(m, color_map_att); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [0,1,2]*1e-3;
ax.XAxis.Visible = 'off';

ylabel('Decoding Strength (norm. µV)'); 
title(tl1, 'Decoding Cu vs Uc', 'FontSize', 8)


%% Fig. C-D. Beh. bias split based on memo maintaence
Behcorr_maint.decoding = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Comb_DistractorProject/CombMemoEnco_cued_600_1100_memarr_September_16_2025_16_00_12.mat');

% DS
x_axis = Behcorr_maint.decoding.comb.lo.wm.bias.x.cu.ds(1,:); 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);

m_lo = mean(Behcorr_maint.decoding.comb.lo.wm.bias.y.cu.ds,1);
m_hi = mean(Behcorr_maint.decoding.comb.hi.wm.bias.y.cu.ds,1);

se_lo = (std(Behcorr_maint.decoding.comb.lo.wm.bias.y.cu.ds,[],1))./sqrt(nEEG);
se_hi = (std(Behcorr_maint.decoding.comb.hi.wm.bias.y.cu.ds,[],1))./sqrt(nEEG);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx, 3);


tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);

hold on; 
ylim([-1, 1]*6)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -20 20 40], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');

% --- p1 (dashed, transparent, open markers) ---
p1_err = errorbar(x_axis(x_idx), k_lo, se_lo(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', color_map(1,:), ...          % open circles
    'MarkerEdgeColor', color_map(1,:), ...
    'Color', color_map(1,:), ...      % errorbar color
    'LineStyle', 'none', ...
    'LineWidth', 1, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p1_line = plot(x_axis, y_dog_lo, '--', ...  % dashed line
    'LineWidth', 2, ...
    'Color', color_map(1,:));   % fully opaque

% --- p2 (solid, opaque, filled markers) ---
p2_err = errorbar(x_axis(x_idx), k_hi, se_hi(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', 'w', ... % filled circles
    'MarkerEdgeColor', color_map(1,:), ...
    'Color', color_map(1,:), ...
    'LineStyle', 'none', ...
    'LineWidth', 0.5, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p2_line = plot(x_axis, y_dog_hi, '-', ...   % solid line
    'LineWidth', 2, ...
    'Color', [color_map(1,:), 0.5]);         % add alpha for transparency

r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YAxis.Visible = 'on';
title('Encoding');
tx = text([15, 15],[-4, -5],{'Weak','Strong'}); 
tx(1).Color = color_map(1,:); tx(2).Color = color_map(2,:);
% ax.XAxis.Visible = 'off'; % remove x-axis
% xlabel('Distractor - Target Orientation (\circ)'); 
% ylabel('Response - Target Orientation (\circ)');



m = [Behcorr_maint.decoding.comb.lo.wm.bias.y.area.uc.ds, Behcorr_maint.decoding.comb.hi.wm.bias.y.area.uc.ds];
doPermutationTest(m(:,1)*[1 0],'larger','right'); doPermutationTest(m(:,2)*[1 0],'larger','right');
doPermutationTest(m,'larger','right');

tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 3;
tl1.Layout.TileSpan = [2 3];
tp = nexttile(tl1, 1);

xlim([0.5 2.5])
ylim([-12, 14]);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;

hold on;
[vp, bp] = violin_plot(m, [color_map(1,:);color_map(1,:)]); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [-1:1:1]*10;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'on';
% tx = text([1,2], [12, 12], {'n.s.', 'p=0.01'}, 'HorizontalAlignment', 'center');
% tx(2).FontSize = 12; 
tx = text([1,2], [-12, -12], {'Weak', 'Strong'}, 'HorizontalAlignment', 'center');
tx(1).Color = color_map(1,:); tx(2).Color = color_map(1,:);

%DO
x_axis = Behcorr_maint.decoding.comb.lo.wm.bias.x.cu.do(1,:); 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);

m_lo = mean(Behcorr_maint.decoding.comb.lo.wm.bias.y.cu.do,1);
m_hi = mean(Behcorr_maint.decoding.comb.hi.wm.bias.y.cu.do,1);

se_lo = (std(Behcorr_maint.decoding.comb.lo.wm.bias.y.cu.do,[],1))./sqrt(nEEG);
se_hi = (std(Behcorr_maint.decoding.comb.hi.wm.bias.y.cu.do,[],1))./sqrt(nEEG);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx, 3);


tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 13;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);

hold on; 
ylim([-1, 1]*6)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -20 20 40], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');

% --- p1 (dashed, transparent, open markers) ---
p1_err = errorbar(x_axis(x_idx), k_lo, se_lo(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', color_map(2,:), ... % filled circles
    'MarkerEdgeColor', color_map(2,:), ...
    'Color', color_map(2,:), ...
    'LineStyle', 'none', ...
    'LineWidth', 0.5, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p1_line = plot(x_axis, y_dog_lo, '--', ...  % dashed line
    'LineWidth', 2, ...
    'Color', color_map(2,:));   % add alpha for transparency

% --- p2 (solid, opaque, filled markers) ---
p2_err = errorbar(x_axis(x_idx), k_hi, se_hi(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', 'w', ...          % open circles
    'MarkerEdgeColor', color_map(2,:), ...
    'Color', color_map(2,:), ...      % errorbar color
    'LineStyle', 'none', ...
    'LineWidth', 1, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p2_line = plot(x_axis, y_dog_hi, '-', ...   % solid line
    'LineWidth', 2, ...
    'Color', [color_map(2,:) 0.5]);         % fully opaque

r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YAxis.Visible = 'on';
title('Encoding');
tx = text([15, 15],[-4, -5],{'Weak','Strong'}); 
tx(1).Color = color_map(2,:); tx(2).Color = color_map(2,:);
% ax.XAxis.Visible = 'off'; % remove x-axis
% xlabel('Distractor - Target Orientation (\circ)'); 
% ylabel('Response - Target Orientation (\circ)');



m = [Behcorr_maint.decoding.comb.lo.wm.bias.y.area.uc.do, Behcorr_maint.decoding.comb.hi.wm.bias.y.area.uc.do];
doPermutationTest(m(:,1)*[1 0],'smaller','left'); doPermutationTest(m(:,2)*[1 0],'smaller','left');
doPermutationTest(m,'smaller','left');

tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 15;
tl1.Layout.TileSpan = [2 3];
tp = nexttile(tl1, 1);

xlim([0.5 2.5])
ylim([-12, 14]);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;

hold on;
[vp, bp] = violin_plot(m, [color_map(2,:);color_map(2,:)]); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [-1:1:1]*10;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'on';
% tx = text([1,2], [12, 12], {'n.s.', 'p=0.01'}, 'HorizontalAlignment', 'center');
% tx(2).FontSize = 12; 
tx = text([1,2], [-12, -12], {'Weak', 'Strong'}, 'HorizontalAlignment', 'center');
tx(1).Color = color_map(2,:); tx(2).Color = color_map(2,:);






%% %%%% memorandum encoding 

%% Fig. XX. Beh split memo maintaence
Behcorr_maint.decoding = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Comb_DistractorProject/CombMemoEnco_cued_100_500_memarr_September_16_2025_19_02_50.mat');

% DS
x_axis = Behcorr_maint.decoding.comb.lo.wm.bias.x.cu.ds(1,:); 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);

m_lo = mean(Behcorr_maint.decoding.comb.lo.wm.bias.y.cu.ds,1);
m_hi = mean(Behcorr_maint.decoding.comb.hi.wm.bias.y.cu.ds,1);

se_lo = (std(Behcorr_maint.decoding.comb.lo.wm.bias.y.cu.ds,[],1))./sqrt(nEEG);
se_hi = (std(Behcorr_maint.decoding.comb.hi.wm.bias.y.cu.ds,[],1))./sqrt(nEEG);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx, 3);


tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);

hold on; 
ylim([-1, 1]*6)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -20 20 40], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');

% --- p1 (dashed, transparent, open markers) ---
p1_err = errorbar(x_axis(x_idx), k_lo, se_lo(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', color_map(1,:), ...          % open circles
    'MarkerEdgeColor', color_map(1,:), ...
    'Color', color_map(1,:), ...      % errorbar color
    'LineStyle', 'none', ...
    'LineWidth', 1, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p1_line = plot(x_axis, y_dog_lo, '-', ...  % dashed line
    'LineWidth', 2, ...
    'Color', color_map(1,:));   % fully opaque

% --- p2 (solid, opaque, filled markers) ---
p2_err = errorbar(x_axis(x_idx), k_hi, se_hi(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', 'w', ... % filled circles
    'MarkerEdgeColor', color_map(1,:), ...
    'Color', color_map(1,:), ...
    'LineStyle', 'none', ...
    'LineWidth', 0.5, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p2_line = plot(x_axis, y_dog_hi, '--', ...   % solid line
    'LineWidth', 2, ...
    'Color', [color_map(1,:), 0.5]);         % add alpha for transparency

r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YAxis.Visible = 'on';
title('Encoding');
tx = text([15, 15],[-4, -5],{'Weak','Strong'}); 
tx(1).Color = color_map(1,:); tx(2).Color = color_map(2,:);
% ax.XAxis.Visible = 'off'; % remove x-axis
% xlabel('Distractor - Target Orientation (\circ)'); 
% ylabel('Response - Target Orientation (\circ)');



m = [Behcorr_maint.decoding.comb.lo.wm.bias.y.area.cu.ds, Behcorr_maint.decoding.comb.hi.wm.bias.y.area.cu.ds];
doPermutationTest(m(:,1)*[1 0],'larger','right'); doPermutationTest(m(:,2)*[1 0],'larger','right');
doPermutationTest(m,'larger','right');

tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 3;
tl1.Layout.TileSpan = [2 3];
tp = nexttile(tl1, 1);

xlim([0.5 2.5])
ylim([-12, 14]);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;

hold on;
[vp, bp] = violin_plot(m, [color_map(1,:);color_map(1,:)]); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [-1:1:1]*10;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'on';
% tx = text([1,2], [12, 12], {'n.s.', 'p=0.01'}, 'HorizontalAlignment', 'center');
% tx(2).FontSize = 12; 
tx = text([1,2], [-12, -12], {'Weak', 'Strong'}, 'HorizontalAlignment', 'center');
tx(1).Color = color_map(1,:); tx(2).Color = color_map(1,:);

%DO
x_axis = Behcorr_maint.decoding.comb.lo.wm.bias.x.cu.do(1,:); 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);

m_lo = mean(Behcorr_maint.decoding.comb.lo.wm.bias.y.cu.do,1);
m_hi = mean(Behcorr_maint.decoding.comb.hi.wm.bias.y.cu.do,1);

se_lo = (std(Behcorr_maint.decoding.comb.lo.wm.bias.y.cu.do,[],1))./sqrt(nEEG);
se_hi = (std(Behcorr_maint.decoding.comb.hi.wm.bias.y.cu.do,[],1))./sqrt(nEEG);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx, 3);


tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 13;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);

hold on; 
ylim([-1, 1]*6)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -20 20 40], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');

% --- p1 (dashed, transparent, open markers) ---
p1_err = errorbar(x_axis(x_idx), k_lo, se_lo(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', color_map(2,:), ... % filled circles
    'MarkerEdgeColor', color_map(2,:), ...
    'Color', color_map(2,:), ...
    'LineStyle', 'none', ...
    'LineWidth', 0.5, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p1_line = plot(x_axis, y_dog_lo, '-', ...  % dashed line
    'LineWidth', 2, ...
    'Color', color_map(2,:));   % add alpha for transparency

% --- p2 (solid, opaque, filled markers) ---
p2_err = errorbar(x_axis(x_idx), k_hi, se_hi(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', 'w', ...          % open circles
    'MarkerEdgeColor', color_map(2,:), ...
    'Color', color_map(2,:), ...      % errorbar color
    'LineStyle', 'none', ...
    'LineWidth', 1, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p2_line = plot(x_axis, y_dog_hi, '--', ...   % solid line
    'LineWidth', 2, ...
    'Color', [color_map(2,:) 0.5]);         % fully opaque

r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YAxis.Visible = 'on';
title('Encoding');
tx = text([15, 15],[-4, -5],{'Weak','Strong'}); 
tx(1).Color = color_map(2,:); tx(2).Color = color_map(2,:);
% ax.XAxis.Visible = 'off'; % remove x-axis
% xlabel('Distractor - Target Orientation (\circ)'); 
% ylabel('Response - Target Orientation (\circ)');



m = [Behcorr_maint.decoding.comb.lo.wm.bias.y.area.cu.do, Behcorr_maint.decoding.comb.hi.wm.bias.y.area.cu.do];
doPermutationTest(m(:,1)*[1 0],'smaller','left'); doPermutationTest(m(:,2)*[1 0],'smaller','left');
doPermutationTest(m,'smaller','left');

tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 15;
tl1.Layout.TileSpan = [2 3];
tp = nexttile(tl1, 1);

xlim([0.5 2.5])
ylim([-12, 14]);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;

hold on;
[vp, bp] = violin_plot(m, [color_map(2,:);color_map(2,:)]); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [-1:1:1]*10;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'on';
% tx = text([1,2], [12, 12], {'n.s.', 'p=0.01'}, 'HorizontalAlignment', 'center');
% tx(2).FontSize = 12; 
tx = text([1,2], [-12, -12], {'Weak', 'Strong'}, 'HorizontalAlignment', 'center');
tx(1).Color = color_map(2,:); tx(2).Color = color_map(2,:);

