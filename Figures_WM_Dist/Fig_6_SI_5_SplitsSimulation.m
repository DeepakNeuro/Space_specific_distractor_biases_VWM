%% Load Simulation data
% Computations based on the model in Scripts/SimulationCodes/Model_main.m

% path settings
addpath(genpath('../Beh_bias_eval/'));

% main model sim
load('../../Data/SimulationsData/Simulations_Comp_March_23_2023_15_53_22.mat');

% sim 2 is memo. maint.  sim1 is dist. encoding
sim_new = load('../../Data/SimulationsData/Simulations_March_23_2023_05_13_56.mat'); % 10_2_1 and 10_2_2

%dist. gating
sim_gat = load('../../Data/SimulationsData/Simulations_VC_HC_delay_0p01_September_28_2025_05_57_52.mat');
% mem. encoding
sim_enco = load('../../Data/SimulationsData/Enco_5x_stim_cu_uc_Simulations_September_28_2025_06_38_13.mat');


%% color maps for various figures
color_map = lines(5);
color_map(6,:) = [0 0 0]; 

color_map_att = brewermap(5, 'PiYG');
color_map_att = color_map_att([1, 5],:); 

color_map_split = color_map(3:4,:); 


%%
plotDefaults; 
tl = tiledlayout(9,9,'TileSpacing','compact', 'Padding', 'compact');


%% Paramters and bias relatioship
x_axis = sim.x; 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
idx = ~isnan(sim.sim.ang_cu_mat(:,1,1)); 

ang_cu = sim.sim.ori_mat(idx,1);
ang_dist = sim.sim.ori_mat(idx,3);

idx_dist_cu = sim.sim.dist_params_mat(idx,2) ==1;

% Sim
ang_resp_cu = circ_grat_mean(squeeze(sim.sim.ang_cu_mat(idx,end,:))')';
ang_diff_cu = circ_dist_grat(ang_resp_cu, ang_cu);

for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim.sim.y_cu_dist_cu(subNo,:), sim.sim.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end



% Sim1
ang_resp_cu = circ_grat_mean(squeeze(sim_new.sim1.ang_cu_mat(idx,end,:))')';
ang_diff_cu = circ_dist_grat(ang_resp_cu, ang_cu);
for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim_new.sim1.y_cu_dist_cu(subNo,:), sim_new.sim1.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end

% Sim2 
ang_resp_cu = circ_grat_mean(squeeze(sim_new.sim2.ang_cu_mat(idx,end,:))')';
ang_diff_cu = circ_dist_grat(ang_resp_cu, ang_cu);

for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim_new.sim2.y_cu_dist_cu(subNo,:), sim_new.sim2.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end

% SimGating
ang_resp_cu = circ_grat_mean(squeeze(sim_gat.sim.ang_cu_mat(idx,end,:))')';
ang_diff_cu = circ_dist_grat(ang_resp_cu, ang_cu);

for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim_gat.sim.y_cu_dist_cu(subNo,:), sim_gat.sim.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end

% Simenco
ang_resp_cu = circ_grat_mean(squeeze(sim_enco.sim_enco.ang_cu_mat(idx,end,:))')';
ang_diff_cu = circ_dist_grat(ang_resp_cu, ang_cu);

for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim_enco.sim_enco.y_cu_dist_cu(subNo,:), sim_enco.sim_enco.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end

%% Distractor time
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [2 2];

tp = nexttile(tl1,1);
x_axis = sim.x; 
m_lo = mean(sim.y_cu_dist_cu_e,1);
m_hi = mean(sim.y_cu_dist_cu_l,1);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

hold on; 
ylim([-1, 1]*3.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), k_lo, m_lo(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
p = errorbar(x_axis(x_idx), k_hi, m_hi(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
plot(x_axis, y_dog_lo, 'LineWidth', 2, 'Color', color_map_split(1,:)); 
plot(x_axis, y_dog_hi, 'LineWidth', 2, 'Color', color_map_split(2,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;
title('Dist. Timing');
tx = text([15, 15],[-2, -3],{'Early','Late'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Stimulus Ori. (°)')
ylabel('Response - Stimulus Ori. (°)')

%% Distractor encoding strength
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 3;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);

m_lo = mean(sim.sim.y_cu_dist_cu,1);
m_hi = mean(sim_new.sim1.y_cu_dist_cu,1);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

hold on; 
ylim([-1, 1]*3.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), k_lo, m_lo(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
p = errorbar(x_axis(x_idx), k_hi, m_hi(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
plot(x_axis, y_dog_lo, 'LineWidth', 2, 'Color', color_map_split(1,:)); 
plot(x_axis, y_dog_hi, 'LineWidth', 2, 'Color', color_map_split(2,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;
% ax.YAxis.Visible = 'off';
title('Dist. Encoding');
tx = text([15, 15],[-2, -3],{'Weak','Strong'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Stimulus Ori. (°)')

% ax.XAxis.Visible = 'off'; % remove x-axis
% xlabel('Distractor - Target Orientation (\circ)'); 
% ylabel('Response - Target Orientation (\circ)'); 

%% Memorandum maintainence
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 7;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);

m_lo = mean(sim_new.sim2.y_cu_dist_cu,1); % stronger maintainence
m_hi = mean(sim.sim.y_cu_dist_cu,1);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

hold on; 
ylim([-1, 1]*3.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), k_lo, m_lo(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
p = errorbar(x_axis(x_idx), k_hi, m_hi(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
plot(x_axis, y_dog_lo, 'LineWidth', 2, 'Color', color_map_split(1,:)); 
plot(x_axis, y_dog_hi, 'LineWidth', 2, 'Color', color_map_split(2,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;
% ax.YAxis.Visible = 'off';
title('Mem. maint.');
tx = text([15, 15],[-2, -3],{'Strong Maint.','Weak Maint.'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);

%% Mem. encoding
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 19;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);


m_lo = mean(sim_enco.sim_enco.y_cu_dist_cu,1);
m_hi = mean(sim.sim.y_cu_dist_cu,1);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

hold on; 
ylim([-1, 1]*3.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), k_lo, m_lo(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
p = errorbar(x_axis(x_idx), k_hi, m_hi(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
plot(x_axis, y_dog_lo, 'LineWidth', 2, 'Color', color_map_split(1,:)); 
plot(x_axis, y_dog_hi, 'LineWidth', 2, 'Color', color_map_split(2,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;
% ax.YAxis.Visible = 'off';
title('Mem. Enco');
tx = text([15, 15],[-2, -3],{'Strong enco.','Weak enco.'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);


%% Dist. Gating
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 5;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);


m_lo = mean(sim_gat.sim.y_cu_dist_cu,1);
m_hi = mean(sim.sim.y_cu_dist_cu,1);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

hold on; 
ylim([-1,1]*3.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), k_lo, m_lo(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
p = errorbar(x_axis(x_idx), k_hi, m_hi(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
plot(x_axis, y_dog_lo, 'LineWidth', 2, 'Color', color_map_split(2,:)); 
plot(x_axis, y_dog_hi, 'LineWidth', 2, 'Color', color_map_split(1,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;
% ax.YAxis.Visible = 'off';
title('Dist. Gating');
tx = text([15, 15],[-2, -3],{'Weak Gat.','Strong Gat.'}); 
tx(1).Color = color_map_split(2,:); tx(2).Color = color_map_split(1,:);


%% %%%%%%%
%{
For "cued" but "distractor opposite" SI 5 splits simulations
%}

%% Paramters and bias relatioship
x_axis = sim.x; 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
idx = ~isnan(sim.sim.ang_cu_mat(:,1,1)); 

ang_cu = sim.sim.ori_mat(idx,1);
ang_dist = sim.sim.ori_mat(idx,3);

idx_dist_cu = sim.sim.dist_params_mat(idx,2) == 2; %dist on uncued side

% Sim
ang_resp_cu = circ_grat_mean(squeeze(sim.sim.ang_cu_mat(idx,end,:))')';
ang_diff_cu = circ_dist_grat(ang_resp_cu, ang_cu);

for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim.sim.y_cu_dist_cu(subNo,:), sim.sim.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end



% Sim1
ang_resp_cu = circ_grat_mean(squeeze(sim_new.sim1.ang_cu_mat(idx,end,:))')';
ang_diff_cu = circ_dist_grat(ang_resp_cu, ang_cu);
for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim_new.sim1.y_cu_dist_cu(subNo,:), sim_new.sim1.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end

% Sim2 
ang_resp_cu = circ_grat_mean(squeeze(sim_new.sim2.ang_cu_mat(idx,end,:))')';
ang_diff_cu = circ_dist_grat(ang_resp_cu, ang_cu);

for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim_new.sim2.y_cu_dist_cu(subNo,:), sim_new.sim2.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end

% SimGating
ang_resp_cu = circ_grat_mean(squeeze(sim_gat.sim.ang_cu_mat(idx,end,:))')';
ang_diff_cu = circ_dist_grat(ang_resp_cu, ang_cu);

for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim_gat.sim.y_cu_dist_cu(subNo,:), sim_gat.sim.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end

% Simenco
ang_resp_cu = circ_grat_mean(squeeze(sim_enco.sim_enco.ang_cu_mat(idx,end,:))')';
ang_diff_cu = circ_dist_grat(ang_resp_cu, ang_cu);

for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim_enco.sim_enco.y_cu_dist_cu(subNo,:), sim_enco.sim_enco.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end

%% Distractor time
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 1+36;
tl1.Layout.TileSpan = [2 2];

tp = nexttile(tl1,1);
x_axis = sim.x; 
m_lo = mean(sim.y_cu_dist_uc_e,1);
m_hi = mean(sim.y_cu_dist_uc_l,1);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

hold on; 
ylim([-1, 1]*3.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), k_lo, m_lo(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
p = errorbar(x_axis(x_idx), k_hi, m_hi(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
plot(x_axis, y_dog_lo, 'LineWidth', 2, 'Color', color_map_split(1,:)); 
plot(x_axis, y_dog_hi, 'LineWidth', 2, 'Color', color_map_split(2,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;
title('Dist. Timing (DO)');
tx = text([15, 15],[-2, -3],{'Early','Late'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Stimulus Ori. (°)')
ylabel('Response - Stimulus Ori. (°)')

%% Distractor encoding strength
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 3+36;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);

m_lo = mean(sim.sim.y_cu_dist_cu,1);
m_hi = mean(sim_new.sim1.y_cu_dist_cu,1);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

hold on; 
ylim([-1, 1]*3.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), k_lo, m_lo(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
p = errorbar(x_axis(x_idx), k_hi, m_hi(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
plot(x_axis, y_dog_lo, 'LineWidth', 2, 'Color', color_map_split(1,:)); 
plot(x_axis, y_dog_hi, 'LineWidth', 2, 'Color', color_map_split(2,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;
% ax.YAxis.Visible = 'off';
title('Dist. Encoding(DO)');
tx = text([15, 15],[-2, -3],{'Weak','Strong'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Stimulus Ori. (°)')

% ax.XAxis.Visible = 'off'; % remove x-axis
% xlabel('Distractor - Target Orientation (\circ)'); 
% ylabel('Response - Target Orientation (\circ)'); 

%% Memorandum maintainence
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 7+36;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);

m_lo = mean(sim_new.sim2.y_cu_dist_cu,1); % stronger maintainence
m_hi = mean(sim.sim.y_cu_dist_cu,1);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

hold on; 
ylim([-1, 1]*3.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), k_lo, m_lo(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
p = errorbar(x_axis(x_idx), k_hi, m_hi(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
plot(x_axis, y_dog_lo, 'LineWidth', 2, 'Color', color_map_split(1,:)); 
plot(x_axis, y_dog_hi, 'LineWidth', 2, 'Color', color_map_split(2,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;
% ax.YAxis.Visible = 'off';
title('Mem. maint.');
tx = text([15, 15],[-2, -3],{'Strong Maint.','Weak Maint.'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);

%% Mem. encoding
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 19+36;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);


m_lo = mean(sim_enco.sim_enco.y_cu_dist_cu,1);
m_hi = mean(sim.sim.y_cu_dist_cu,1);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

hold on; 
ylim([-1, 1]*3.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), k_lo, m_lo(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
p = errorbar(x_axis(x_idx), k_hi, m_hi(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
plot(x_axis, y_dog_lo, 'LineWidth', 2, 'Color', color_map_split(1,:)); 
plot(x_axis, y_dog_hi, 'LineWidth', 2, 'Color', color_map_split(2,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;
% ax.YAxis.Visible = 'off';
title('Mem. Enco(DO)');
tx = text([15, 15],[-2, -3],{'Strong enco.','Weak enco.'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);


%% Dist. Gating
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 5+36;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);


m_lo = mean(sim_gat.sim.y_cu_dist_cu,1);
m_hi = mean(sim.sim.y_cu_dist_cu,1);

k_lo = warp_mov_mean(m_lo(x_idx),3);
k_hi = warp_mov_mean(m_hi(x_idx),3);

y_dog_lo = fit_curves(x_axis(x_idx), k_lo, x_idx);
y_dog_hi = fit_curves(x_axis(x_idx), k_hi, x_idx);

hold on; 
ylim([-1,1]*3.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), k_lo, m_lo(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
p = errorbar(x_axis(x_idx), k_hi, m_hi(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4);
plot(x_axis, y_dog_lo, 'LineWidth', 2, 'Color', color_map_split(2,:)); 
plot(x_axis, y_dog_hi, 'LineWidth', 2, 'Color', color_map_split(1,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;
% ax.YAxis.Visible = 'off';
title('Dist. Gating(DO)');
tx = text([15, 15],[-2, -3],{'Weak Gat.','Strong Gat.'}); 
tx(1).Color = color_map_split(2,:); tx(2).Color = color_map_split(1,:);

 