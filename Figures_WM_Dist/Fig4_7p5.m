%{
Code for reproducing Fig. 4 of the paper. 
%}
Figures_main

%% Load Data
Behcorr.time = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Comb_DistractorProject/Comb_March_24_2023_09_00_55.mat');
Behcorr.gating = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Comb_DistractorProject/Comb_May_21_2024_15_27_19.mat');
Behcorr.decoding = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Comb_DistractorProject/Comb_May_21_2024_16_17_19.mat');

%% Fig. 4 A distractor evoked ERP
eeg.erp = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/EEG/Decoding/EEGAnalysis_ERPsDist_May_17_2024_13_45_32.mat');

e_c_ad = eeg.erp.eeg.erp.distractorOnset.e_c_ad;
e_c_ud = eeg.erp.eeg.erp.distractorOnset.e_c_ud;
e_i_ad = eeg.erp.eeg.erp.distractorOnset.e_i_ad;
e_i_ud = eeg.erp.eeg.erp.distractorOnset.e_i_ud;
time_erp = eeg.erp.eeg.erp.distractorOnset.time(1,:); 

m_c_cu = mean(e_c_ad,1);
m_c_uc = mean(e_c_ud,1);
m_i_cu = mean(e_i_ad,1);
m_i_uc = mean(e_i_ud,1);

se_c_cu = std(e_c_ad,[],1)./sqrt(nEEG);
se_c_uc = std(e_c_ud,[],1)./sqrt(nEEG);
se_i_cu = std(e_i_ad,[],1)./sqrt(nEEG);
se_i_uc = std(e_i_ud,[],1)./sqrt(nEEG);

m_c = [m_c_cu; m_c_uc]';
se_c = [se_c_cu; se_c_uc]';
m_i = [m_i_cu; m_i_uc]';
se_i = [se_i_cu; se_i_uc]';

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact'); 
tl1.Layout.Tile = 19;
tl1.Layout.TileSpan = [2 3];

tp = nexttile(tl1,1);

hold on; 
ylim([-1.5, 1.1])
xlim([-0.05, 0.4])
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0.15); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0.2); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0.35); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = boundedline(time_erp, m_c , se_c, 'cmap', color_map_att);
ax = gca; 
ax.Box = 'off';
ax.XTick = [0:0.1:0.4];
ax.YTick = [-1:1:1];
ylabel('Contra ERP (µV)')

% ERP Quantification
%% Find the peaks in the ERP (200- 270 & 270 - 350)
toi_p2 = (time_erp>0.2 & time_erp<0.27);
toi_p3a = (time_erp>0.27 & time_erp<0.35);
for str = ["ad","ud"]
    eval(sprintf('[~,peak_idx_toi_p2] = max(e_c_%s(:, toi_p2),[],2);', str));
    eval(sprintf('[~,peak_idx_toi_p3a] = max(e_c_%s(:, toi_p3a),[],2);', str));
    peak_idx_p2 = find(toi_p2,1) + peak_idx_toi_p2-1;
    peak_idx_p3a = find(toi_p3a,1) + peak_idx_toi_p3a-1;
    
    for sub = 1:nEEG
        eval(sprintf('p2_erp_val_c_%s(sub) = e_c_%s(sub, peak_idx_p2(sub));',str, str));
        eval(sprintf('p3a_erp_val_c_%s(sub) = e_c_%s(sub, peak_idx_p3a(sub));', str, str));
    end
end
p2_p3a_comb_cu = (p2_erp_val_c_ad + p3a_erp_val_c_ad)./2;
p2_p3a_comb_uc = (p2_erp_val_c_ud + p3a_erp_val_c_ud)./2;
%%
t_idx_n1 = (time_erp >=0.15)&(time_erp <=0.2); % 0.15 to 0.2 s from dist onset    
t_idx_p2 = (time_erp >=0.2)&(time_erp <= 0.35); % 0.2 to 0.35 s from dist onset 

erp_cu_contra_n1 = mean(e_c_ad(:,t_idx_n1),2);
erp_uc_contra_n1 = mean(e_c_ud(:,t_idx_n1),2);

erp_cu_ipsi_n1 = mean(e_i_ad(:,t_idx_n1),2);
erp_uc_ipsi_n1 = mean(e_i_ud(:,t_idx_n1),2);

erp_cu_contra_p2 = p2_p3a_comb_cu.';
erp_uc_contra_p2 = p2_p3a_comb_uc.';

erp_cu_ipsi_p2 = mean(e_i_ad(:,t_idx_p2),2);
erp_uc_ipsi_p2 = mean(e_i_ud(:,t_idx_p2),2);

m_n1 = [erp_cu_contra_n1, erp_uc_contra_n1]; %, erp_cu_ipsi_n1, erp_uc_ipsi_n1];
m_p2 = [erp_cu_contra_p2, erp_uc_contra_p2]; %, erp_cu_ipsi_p2, erp_uc_ipsi_p2];

doPermutationTest([erp_cu_contra_n1, erp_uc_contra_n1],'both','both')
doPermutationTest([erp_cu_contra_p2, erp_uc_contra_p2],'both','both')

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact'); 
tl1.Layout.Tile = 22;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1,1);
xlim([0.5 2.5])
ylim([-4 2]);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
hold on;
[vp, bp] = violin_plot(m_n1, color_map_att); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [-3,0];
ax.XAxis.Visible = 'off';
title(tl1, 'ERP Cu vs Uc (N2)', 'FontSize', 8)

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact'); 
tl1.Layout.Tile = 34;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1,1);
xlim([0.5 2.5])
ylim([-2 3.2]);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
hold on;
[vp, bp] = violin_plot(m_p2, color_map_att); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [0, 2];
ax.XAxis.Visible = 'off';

title(tl1, 'ERP Cu vs Uc (p2/p3a)', 'FontSize', 8)



%% Fig. 4B & corresponding SI. Median split Distractor time
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);
x_axis = Behcorr.decoding.comb.lo.wm.bias.x.cu.ds(1,:); 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
m_lo = nanmean(Behcorr.time.comb.lo.wm.bias.y.cu.ds,1);
m_hi = nanmean(Behcorr.time.comb.hi.wm.bias.y.cu.ds,1);

se_lo = (nanstd(Behcorr.time.comb.lo.wm.bias.y.cu.ds,[],1))./sqrt(nSubjects);
se_hi = (nanstd(Behcorr.time.comb.hi.wm.bias.y.cu.ds,[],1))./sqrt(nSubjects);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);


y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx); %y_dog_lo = -y_dog_lo; 
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

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
    'LineWidth', 1, ...
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
ax.YTick = [-3:1:3]*2;
title('Time');
tx = text([15, 15],[-4, -5],{'Early','Late'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Actual Ori. (°)')
ylabel('Response - Actual Ori. (°)')


m = [Behcorr.time.comb.lo.wm.bias.y.area.cu.ds, Behcorr.time.comb.hi.wm.bias.y.area.cu.ds];
% doPermutationTest(m(:,1)*[1 0],'smaller'); doPermutationTest(m(:,2)*[1 0],'smaller');
% doPermutationTest(m,'larger');

% violin plot
tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 13;
tl1.Layout.TileSpan = [2 3];
tp = nexttile(tl1, 1);

xlim([0.5 2.5])
ylim([-12, 14]);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
hold on;
[vp, bp] = violin_plot(m,  [color_map(1,:); color_map(1,:)]); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [-1:1:1]*10;
ax.XAxis.Visible = 'off';
ylabel('Bias (°)'); 
% tx = text([1,2], [12, 12], {'n.s.', 'p=0.001'}, 'HorizontalAlignment', 'center');
% tx(2).FontSize = 12; 
tx = text([1,2], [-12, -12], {'Early', 'Late'}, 'HorizontalAlignment', 'center');
tx(1).Color = color_map(1,:); tx(2).Color = color_map(1,:);


%% Fig. 4C. & corresponding SI --  Median split Distractor gating 
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 5;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);

m_lo = nanmean(Behcorr.gating.comb.hi.wm.bias.y.cu.ds,1); % High P2 is low gating
m_hi = nanmean(Behcorr.gating.comb.lo.wm.bias.y.cu.ds,1); % Low P2 is high gating

se_lo = (nanstd(Behcorr.gating.comb.hi.wm.bias.y.cu.ds,[],1))./sqrt(nEEG);
se_hi = (nanstd(Behcorr.gating.comb.lo.wm.bias.y.cu.ds,[],1))./sqrt(nEEG);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

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
    'LineWidth', 1, ...
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
ax.Ytick = [-6:2:6];
ax.YAxis.Visible = 'on';
title('Gating');
tx = text([15, 15],[-4, -5],{'Strong P2/P3a','Weak P2/P3a'}); 
tx(1).Color = color_map(1,:); tx(2).Color = color_map(1,:);
xlabel('Distractor - Actual Ori. (°)')



non_nan_hi = ~isnan(Behcorr.gating.comb.hi.wm.bias.y.area.cu.do);
non_nan_lo = ~isnan(Behcorr.gating.comb.lo.wm.bias.y.area.cu.do);
m = [Behcorr.gating.comb.hi.wm.bias.y.area.cu.do(non_nan_hi), Behcorr.gating.comb.lo.wm.bias.y.area.cu.do(non_nan_lo)];
% doPermutationTest(m(:,1)*[1 0],'smaller'); doPermutationTest(m(:,2)*[1 0],'smaller');
% doPermutationTest(m,'larger');

tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 17;
tl1.Layout.TileSpan = [2 3];
tp = nexttile(tl1, 1);

xlim([0.5 2.5])
ylim([-12, 14]);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
hold on; 
[vp, bp] = violin_plot(m, color_map_split); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [-1:1:1]*10;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'on';
% tx = text([1,2], [12, 12], {'n.s.', 'p=0.002'}, 'HorizontalAlignment', 'center');
% tx(2).FontSize = 12; 
tx = text([1,2], [-12, -12], {'Strong P2/P3a', 'Weak P2/P3a'}, 'HorizontalAlignment', 'center');
tx(1).Color = color_map(1,:); tx(2).Color = color_map(1,:);
