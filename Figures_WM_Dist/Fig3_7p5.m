%{
Code for reproducing Fig. 3 of the paper. 
%}

Figures_main
%% Load distractor decoding output
eeg.dec.dist = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/EEG_posterior_all/EEGAnalysis_Distractor_Decoding_July_06_2024_14_09_56.mat');

%x axis for tuning curves
step_size = 180/16;
angspace_ang2 = (90:-step_size:-90) - step_size/2;
angspace_ang2(end) = [];

%% Fig. 3A. Distractor tuning curve
time = eeg.dec.dist.eeg.decodingTrC.DistractorOnset.time(1,:); 
dis_mem_l = nanmean(eeg.dec.dist.eeg.decodingTrC.DistractorOnset.l.dists_mem, 4);
dis_mem_r = nanmean(eeg.dec.dist.eeg.decodingTrC.DistractorOnset.r.dists_mem, 4);

dis_mem = squeeze(mean((dis_mem_l + dis_mem_r)/2,1));

% figure
color_map_dec = flip(brewermap(250, 'RdBu')); %linspecer(250);

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact'); 
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [2 3];

tp = nexttile(tl1, 1);

imagesc(time, angspace_ang2, dis_mem); box off; 
colormap(color_map_dec); cb = colorbar; caxis([-6 6]*1e-3);
xlabel('Time from distractor onset (s)')
ylabel('Angle Difference (°)')
hold on; 
xlim([-0.1, 0.57])
ylim([-1, 1]*90)
ax = gca; 
ax.Box = 'off';
ax.XTick = [0:0.25:0.5];
ax.YTick = [-90:45:90];


%% Fig. 3B. Distractor Cued vs Uncued decoding 
% Distractor 
time_dist = eeg.dec.dist.eeg.decodingTrC.DistractorOnset.time(1,:); 
cos_mem_l_al = nanmean(eeg.dec.dist.eeg.decodingTrC.DistractorOnset.l.cos_mem(:,:,idxL), 3);
cos_mem_r_al = nanmean(eeg.dec.dist.eeg.decodingTrC.DistractorOnset.r.cos_mem(:,:,idxL), 3);
cos_mem_l_ar = nanmean(eeg.dec.dist.eeg.decodingTrC.DistractorOnset.l.cos_mem(:,:,idxR), 3);
cos_mem_r_ar = nanmean(eeg.dec.dist.eeg.decodingTrC.DistractorOnset.r.cos_mem(:,:,idxR), 3);

cos_mem_avg_cu = (cos_mem_l_al + cos_mem_r_ar)/2;
cos_mem_avg_uc = (cos_mem_l_ar + cos_mem_r_al)/2;
cos_mem_dist_cu = smoothdata(mean(cos_mem_avg_cu,1),'movmean', 10);
cos_mem_dist_uc = smoothdata(mean(cos_mem_avg_uc,1),'movmean', 10);
cos_mem_dist_se_cu = std(cos_mem_avg_cu,[],1)./sqrt(nEEG);
cos_mem_dist_se_uc = std(cos_mem_avg_uc,[],1)./sqrt(nEEG);
m = [cos_mem_dist_cu; cos_mem_dist_uc]';
se = [cos_mem_dist_se_cu; cos_mem_dist_se_uc]';

[pVal_line_dist_cu, pVal_dist_cu] = findSignificanceValueTs_ClusterPermutest(cos_mem_avg_cu, cos_mem_avg_cu*0, true);
[pVal_line_dist_uc, pVal_dist_uc] = findSignificanceValueTs_ClusterPermutest(cos_mem_avg_uc, cos_mem_avg_uc*0, true);

%figure
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact'); 
tl1.Layout.Tile = 13;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1,1);

hold on; 
ylim([-4, 15]*1e-4)
xlim([-0.1, 0.57])
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = boundedline(time_dist, m , se, 'cmap', color_map_att);
plot(time_dist, pVal_line_dist_cu*(-2*1e-4), 'Color', color_map(6,:), 'LineWidth', 1);
plot(time_dist, pVal_line_dist_uc*(-3*1e-4), 'Color', color_map(6,:), 'LineWidth', 1);

ax = gca; 
ax.Box = 'off';
ax.XTick = [0:0.25:0.5];
ax.YTick = [0:6:12]*1e-4;
xlabel('Time (s) from distractor onset');

% comparision and violin plot at peak Dist. decoding
t_idx = knnsearch(time_dist', 0.25); % 166
m = [cos_mem_avg_cu(:,t_idx), cos_mem_avg_uc(:,t_idx)];
doPermutationTest(m*1e4);

% figure
tl1 = tiledlayout(tl, 2,2, 'TileSpacing','compact', 'Padding', 'compact'); 
tl1.Layout.Tile = 5;
tl1.Layout.TileSpan = [2 2];

tp = nexttile(tl1,1);

xlim([0.5 2.5])
ylim([-3 6]*1e-3);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
hold on;
[vp, bp] = violin_plot(m, color_map_att); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [0,3]*1e-3;
ax.XAxis.Visible = 'off';

ylabel('Decoding Strength (norm. µV)'); 
title(tl1, 'Decoding Cu vs Uc', 'FontSize', 8)


%% Fig 3 C-D Distractor Encoding based Median Spits
Behcorr.decoding = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Comb_DistractorProject/Comb_May_21_2024_16_17_19.mat');


tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 25;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);
x_axis = Behcorr.decoding.comb.lo.wm.bias.x.cu.ds(1,:);
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
m_lo = mean(Behcorr.decoding.comb.lo.wm.bias.y.cu.ds,1);
m_hi = mean(Behcorr.decoding.comb.hi.wm.bias.y.cu.ds,1);

se_lo = (std(Behcorr.decoding.comb.lo.wm.bias.y.cu.ds,[],1))./sqrt(nEEG);
se_hi = (std(Behcorr.decoding.comb.hi.wm.bias.y.cu.ds,[],1))./sqrt(nEEG);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx, 3);

hold on; 
ylim([-1, 1]*6)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -20 20 40], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');

% --- p1 (dashed, transparent, open markers) ---
p1_err = errorbar(x_axis(x_idx), k_lo, se_lo(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', 'none', ...          % open circles
    'MarkerEdgeColor', color_map(1,:), ...
    'Color', color_map(1,:), ...      % errorbar color
    'LineStyle', 'none', ...
    'LineWidth', 0.5, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p1_line = plot(x_axis, y_dog_lo, '--', ...  % dashed line
    'LineWidth', 2, ...
    'Color', [color_map(1,:) 0.5]);   % add alpha for transparency

% --- p2 (solid, opaque, filled markers) ---
p2_err = errorbar(x_axis(x_idx), k_hi, se_hi(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', color_map(1,:), ... % filled circles
    'MarkerEdgeColor', color_map(1,:), ...
    'Color', color_map(1,:), ...
    'LineStyle', 'none', ...
    'LineWidth', 0.5, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p2_line = plot(x_axis, y_dog_hi, '-', ...   % solid line
    'LineWidth', 2, ...
    'Color', color_map(1,:));         % fully opaque

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

m = [Behcorr.decoding.comb.lo.wm.bias.y.area.cu.ds, Behcorr.decoding.comb.hi.wm.bias.y.area.cu.ds];
doPermutationTest(m(:,1)*[1 0],'larger', 'right'); doPermutationTest(m(:,2)*[1 0],'larger','right');
doPermutationTest(m,'smaller','left');

tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 27;
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




% DO
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 37;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);
x_axis = Behcorr.decoding.comb.lo.wm.bias.x.cu.ds(1,:); 
m_lo = mean(Behcorr.decoding.comb.lo.wm.bias.y.cu.do,1);
m_hi = mean(Behcorr.decoding.comb.hi.wm.bias.y.cu.do,1);

se_lo = (std(Behcorr.decoding.comb.lo.wm.bias.y.cu.do,[],1))./sqrt(nEEG);
se_hi = (std(Behcorr.decoding.comb.hi.wm.bias.y.cu.do,[],1))./sqrt(nEEG);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx, 3);

hold on; 
ylim([-1, 1]*6)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -20 20 40], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');

% --- p1 (dashed, transparent, open markers) ---
p1_err = errorbar(x_axis(x_idx), k_lo, se_lo(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', 'none', ...          % open circles
    'MarkerEdgeColor', color_map(2,:), ...
    'Color', color_map(2,:), ...      % errorbar color
    'LineStyle', 'none', ...
    'LineWidth', 0.5, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p1_line = plot(x_axis, y_dog_lo, '--', ...  % dashed line
    'LineWidth', 2, ...
    'Color', [color_map(2,:) 0.5]);   % add alpha for transparency

% --- p2 (solid, opaque, filled markers) ---
p2_err = errorbar(x_axis(x_idx), k_hi, se_hi(x_idx), ...
    'o', ...                                % circle markers
    'MarkerFaceColor', color_map(2,:), ... % filled circles
    'MarkerEdgeColor', color_map(2,:), ...
    'Color', color_map(2,:), ...
    'LineStyle', 'none', ...
    'LineWidth', 0.5, ...
    'CapSize', 0, ...
    'MarkerSize', 4);

p2_line = plot(x_axis, y_dog_hi, '-', ...   % solid line
    'LineWidth', 2, ...
    'Color', color_map(2,:));         % fully opaque

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

m = [Behcorr.decoding.comb.lo.wm.bias.y.area.cu.do, Behcorr.decoding.comb.hi.wm.bias.y.area.cu.do];
doPermutationTest(m(:,1)*[1 0],'smaller', 'left'); doPermutationTest(m(:,2)*[1 0],'smaller','left');
doPermutationTest(m,'larger','right');

tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 39;
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
