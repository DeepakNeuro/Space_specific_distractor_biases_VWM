
%{
Code for reporducing figure 6 and SI 5;

1. Code for reproducing 6A & SI 5A of the paper.
2. Code for reproducing 6B-F & SI 5B-F of the paper.
%}
% Figures_main

sim = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Simulations_DistractorProject/Simulations_Comp_March_23_2023_15_53_22.mat');

N = 512; % Number of neurons 
theta_neuron = rad2deg([0:N-1]/N*pi - pi/2);
t_sim = [0:2:3050-2]/1e3 - 0.5; % Time locked to stimulus onset
t_sim_idx = t_sim > 0.1; 

nSim = 10; % 10 groups of 1000 simulations = 10,000 simulations

t_idx = knnsearch(sim.time_rel', 100);
sim.y_po_dist_ds = (sim.y_uc_dist_uc + sim.y_cu_dist_cu)./2;
sim.y_po_dist_do = (sim.y_uc_dist_cu + sim.y_cu_dist_uc)./2;

%% Fig. 6A-B. Model's Behavioral Bias
x_axis = sim.x; 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
m_ds = mean(sim.y_po_dist_ds,1);
m_do = mean(sim.y_po_dist_do,1);

y_dog_ds = fit_curves(x_axis(x_idx), m_ds(x_idx), x_idx);
y_dog_do = fit_curves(x_axis(x_idx), m_do(x_idx), x_idx);

% figure
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);
hold on; 
ylim([-1, 1]*2.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), m_ds(x_idx), m_ds(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0);
p.MarkerSize = 4; 
plot(x_axis, y_dog_ds, 'LineWidth', 2, 'Color', color_map(1,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*2;

xlabel('Distractor - Target Orientation (\circ)'); 
ylabel('Response - Target Orientation (\circ)'); 


tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 13;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);

hold on; 
ylim([-1, 1]*2.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), m_do(x_idx), m_do(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0);
p.MarkerSize = 4; 
plot(x_axis, y_dog_do, 'LineWidth', 2, 'Color', color_map(2,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*2;

xlabel(tl1, 'Distractor - Stimulus Ori. (°)', 'FontSize', 8); 
ylabel(tl1, 'Response - Stimulus Ori. (°)', 'FontSize', 8); 
title(tl1, 'Behavioral Bias', 'FontSize', 8)

%% Fig. 6 C-D. Model's Neural Bias
sim.y_val_po_dist_ds = (sim.y_val_cu_dist_cu+sim.y_val_uc_dist_uc)./2;
sim.y_val_po_dist_do = (sim.y_val_cu_dist_uc+sim.y_val_uc_dist_cu)./2;

x_axis = sim.x; 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
m_ds = mean(squeeze(sim.y_val_po_dist_ds(:,t_idx,:)),1);
m_do = mean(squeeze(sim.y_val_po_dist_do(:,t_idx,:)),1);

y_dog_ds = fit_curves(x_axis(x_idx), m_ds(x_idx), x_idx);
y_dog_do = fit_curves(x_axis(x_idx), m_do(x_idx), x_idx);


% figure
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 3;
tl1.Layout.TileSpan = [2 2];

tp = nexttile(tl1, 1);

hold on; 
ylim([-1, 1]*25)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), m_ds(x_idx), m_ds(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0);
p.MarkerSize = 4; 
plot(x_axis, y_dog_ds, 'LineWidth', 2, 'Color', color_map(1,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*20;
title('Pooled');
xlabel('Distractor - Target Orientation (\circ)'); 
ylabel('Response - Target Orientation (\circ)'); 


tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 15;
tl1.Layout.TileSpan = [2 2];

tp = nexttile(tl1, 1);

hold on; 
ylim([-1, 1]*5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), m_do(x_idx), m_do(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0);
p.MarkerSize = 4; 
plot(x_axis, y_dog_do, 'LineWidth', 2, 'Color', color_map(2,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*4;






%% Load Simulation data
% Computations based on 10_2
addpath '/media/hdd/Sanchit/Exogenous_Project/Data Analysis Codes'

% normal simulation
sim = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Simulations_DistractorProject/Simulations_Comp_March_23_2023_15_53_22.mat');

% memo. maint. in sim2 of this.
sim_new = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Simulations_DistractorProject/Simulations_March_23_2023_05_13_56.mat'); % 10_2_1 and 10_2_2

sim_gat = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Simulations_DistractorProject/Simulations_VC_HC_delay_0p01_September_28_2025_05_57_52.mat');
sim_enco = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Simulations_DistractorProject/Enco_5x_stim_cu_uc_Simulations_September_28_2025_06_38_13.mat');



%% color maps for various figures
color_map = lines(5);
color_map(6,:) = [0 0 0]; 

color_map_att = brewermap(5, 'PiYG');
color_map_att = color_map_att([1, 5],:); 

color_map_split = color_map(3:4,:); 


%%
plotDefaults; 
tl = tiledlayout(9,9,'TileSpacing','compact', 'Padding', 'compact');


%% Behavior bias model
N = 512; % Number of neurons 
theta_neuron = rad2deg([0:N-1]/N*pi - pi/2);
t_sim = [0:2:3050-2]/1e3 - 0.5; % Time locked to stimulus onset
t_sim_idx = t_sim > 0.1; 

nSim = 10; % 10 groups of 1000 simulations = 10,000 simulations

t_idx = knnsearch(sim.time_rel', 100);
sim.y_po_dist_ds = (sim.y_uc_dist_uc + sim.y_cu_dist_cu)./2;
sim.y_po_dist_do = (sim.y_uc_dist_cu + sim.y_cu_dist_uc)./2;

%% Fig. 6A-B. Model's Behavioral Bias
x_axis = sim.x; 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
m_ds = mean(sim.y_po_dist_ds,1);
m_do = mean(sim.y_po_dist_do,1);

y_dog_ds = fit_curves(x_axis(x_idx), m_ds(x_idx), x_idx);
y_dog_do = fit_curves(x_axis(x_idx), m_do(x_idx), x_idx);

% figure
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [3, 3];

tp = nexttile(tl1, 1);
hold on; 
ylim([-1, 1]*2.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), m_ds(x_idx), m_ds(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0);
p.MarkerSize = 4; 
plot(x_axis, y_dog_ds, 'LineWidth', 2, 'Color', color_map(1,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*2;

xlabel('Distractor - Target Orientation (\circ)'); 
ylabel('Response - Target Orientation (\circ)'); 


tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 46;
tl1.Layout.TileSpan = [3, 3];

tp = nexttile(tl1, 1);

hold on; 
ylim([-1, 1]*2.5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -30 20 60], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
p = errorbar(x_axis(x_idx), m_do(x_idx), m_do(x_idx)*0, 'MarkerFaceColor', [1 1 1], 'Color', color_map(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0);
p.MarkerSize = 4; 
plot(x_axis, y_dog_do, 'LineWidth', 2, 'Color', color_map(2,:)); 
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*2;

xlabel(tl1, 'Distractor - Stimulus Ori. (°)', 'FontSize', 8); 
ylabel(tl1, 'Response - Stimulus Ori. (°)', 'FontSize', 8); 
title(tl1, 'Behavioral Bias', 'FontSize', 8)


%% Paramters and bias relatioship
x_axis = sim.x; 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
idx = ~isnan(sim.sim.ang_cu_mat(:,1,1)); 

ang_cu = sim.sim.ori_mat(idx,1);
ang_dist = sim.sim.ori_mat(idx,3);

idx_dist_cu = sim.sim.dist_params_mat(idx,2) ==1;

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
tl1.Layout.Tile = 49;
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
title('Timing');
tx = text([15, 15],[-2, -3],{'Early','Late'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Stimulus Ori. (°)')
ylabel('Response - Stimulus Ori. (°)')

%% Distractor encoding strength
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 8;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);

m_lo = mean(sim.y_cu_dist_cu,1);
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
title('Encoding');
tx = text([15, 15],[-2, -3],{'Weak','Strong'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Stimulus Ori. (°)')

% ax.XAxis.Visible = 'off'; % remove x-axis
% xlabel('Distractor - Target Orientation (\circ)'); 
% ylabel('Response - Target Orientation (\circ)'); 

%% Memorandum maintainence
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 4;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);

m_lo = mean(sim_new.sim2.y_cu_dist_cu,1); % stronger maintainence
m_hi = mean(sim.y_cu_dist_cu,1);

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
ax.XTick = [-90:45:90];sim.y_cu_dist_cu
ax.YTick = [-1:1:1]*3;
% ax.YAxis.Visible = 'off';
title('Mem. maint.');
tx = text([15, 15],[-2, -3],{'Strong Maint.','Weak Maint.'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);

%% Mem. encoding
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 6;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);


m_lo = mean(sim_enco.sim_enco.y_cu_dist_cu,1);
m_hi = mean(sim.y_cu_dist_cu,1);

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
tl1.Layout.Tile = 51;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);


m_lo = mean(sim_gat.sim.y_cu_dist_cu,1);
m_hi = mean(sim.y_cu_dist_cu,1);

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
title('Gating');
tx = text([15, 15],[-2, -3],{'Weak Gat.','Strong Gat.'}); 
tx(1).Color = color_map_split(2,:); tx(2).Color = color_map_split(1,:);


%% %%%%%%%
%{
Code for reproducing 6E-G of the paper. (DO)
%}
%% Paramters and bias relatioship
x_axis = sim.x; 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
idx = ~isnan(sim.sim.ang_cu_mat(:,1,1)); 

ang_cu = sim.sim.ori_mat(idx,1);
ang_dist = sim.sim.ori_mat(idx,3);

idx_dist_cu = sim.sim.dist_params_mat(idx,2) == 2; % dist on the uncued side

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
    [~, ~, ~, sim_gat.sim.y_cu_dist_cu(subNo,:), sim_gat.sim2.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
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
tl1.Layout.Tile = 67;
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
title('Timing');
tx = text([15, 15],[-2, -3],{'Early','Late'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Stimulus Ori. (°)')
ylabel('Response - Stimulus Ori. (°)')

%% Distractor encoding strength
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 35;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);

m_lo = mean(sim.y_cu_dist_uc,1);
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
title('Encoding');
tx = text([15, 15],[-2, -3],{'Weak','Strong'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Stimulus Ori. (°)')

% ax.XAxis.Visible = 'off'; % remove x-axis
% xlabel('Distractor - Target Orientation (\circ)'); 
% ylabel('Response - Target Orientation (\circ)'); 

%% Memorandum maintainence
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 31;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);


m_lo = mean(sim_new.sim2.y_cu_dist_cu,1); % strong maint
m_hi = mean(sim.y_cu_dist_uc,1);

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
ax.XTick = [-90:45:90];sim.y_cu_dist_cu  
ax.YTick = [-1:1:1]*3;
% ax.YAxis.Visible = 'off';
title('Mem. maint.');
tx = text([15, 15],[-2, -3],{'Strong Maint.','Weak Maint.'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);

%% Mem. encoding
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 33;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);


m_lo = mean(sim_enco.sim_enco.y_cu_dist_cu,1);
m_hi = mean(sim.y_cu_dist_uc,1);

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
ax.XTick = [-90:45:90];sim.y_cu_dist_cu  
ax.YTick = [-1:1:1]*3;
% ax.YAxis.Visible = 'off';
title('Mem. Enco');
tx = text([15, 15],[-2, -3],{'Strong Gat.','Weak Gat.'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);




%% Dist. Gating
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 69;
tl1.Layout.TileSpan = [2 2];
tp = nexttile(tl1, 1);


m_lo = mean(sim_gat.sim.y_cu_dist_cu,1);
m_hi = mean(sim.y_cu_dist_uc,1);

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
title('Gating');
tx = text([15, 15],[-2, -3],{'Weak Gat.','Strong Gat.'}); 
tx(1).Color = color_map_split(2,:); tx(2).Color = color_map_split(1,:);


 