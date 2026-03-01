
%{
Code for reporducing figure 6A and SI 5 A;
%}
Figures_main
load('../../Data/SimulationsData/Simulations_Comp_March_23_2023_15_53_22.mat');

N = 512; % Number of neurons 
theta_neuron = rad2deg([0:N-1]/N*pi - pi/2);
t_sim = [0:2:3050-2]/1e3 - 0.5; % Time locked to stimulus onset
t_sim_idx = t_sim > 0.1; 

nSim = 10; % 10 groups of 1000 simulations = 10,000 simulations

t_idx = knnsearch(sim.time_rel', 100);

x_axis = sim.x; 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
idx = ~isnan(sim.sim.ang_cu_mat(:,1,1)); 

ang_cu = sim.sim.ori_mat(idx,1);
ang_uc = sim.sim.ori_mat(idx,2);
ang_dist = sim.sim.ori_mat(idx,3);

idx_dist_cu = sim.sim.dist_params_mat(idx,2) ==1;
idx_dist_uc = sim.sim.dist_params_mat(idx,2) ==2;

% Sim
ang_resp_cu = circ_grat_mean(squeeze(sim.sim.ang_cu_mat(idx,end,:))')';
ang_diff_cu = circ_dist_grat(ang_resp_cu, ang_cu);

for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim.sim.y_cu_dist_cu(subNo,:), sim.sim.area_m_cu_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_cu & tr_idx), ang_cu(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end

for subNo = 1:10
    tr_idx = false(size(ang_cu)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim.sim.y_cu_dist_uc(subNo,:), sim.sim.area_m_cu_dist_uc(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_cu(idx_dist_uc & tr_idx), ang_cu(idx_dist_uc & tr_idx), ang_dist(idx_dist_uc & tr_idx));
end


ang_resp_uc = circ_grat_mean(squeeze(sim.sim.ang_uc_mat(idx,end,:))')';
ang_diff_uc = circ_dist_grat(ang_resp_uc, ang_uc);
for subNo = 1:10
    tr_idx = false(size(ang_uc)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim.sim.y_uc_dist_uc(subNo,:), sim.sim.area_m_uc_dist_uc(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_uc(idx_dist_uc & tr_idx), ang_uc(idx_dist_uc & tr_idx), ang_dist(idx_dist_uc & tr_idx));
end

for subNo = 1:10
    tr_idx = false(size(ang_uc)); 
    tr_idx([1:1000]+(subNo-1)*1000) = true; 
    [~, ~, ~, sim.sim.y_uc_dist_cu(subNo,:), sim.sim.area_m_uc_dist_cu(subNo,:)] = WMA_biasAnalysis_bias(ang_diff_uc(idx_dist_cu & tr_idx), ang_uc(idx_dist_cu & tr_idx), ang_dist(idx_dist_cu & tr_idx));
end

sim.y_po_dist_ds = (sim.sim.y_cu_dist_cu + sim.sim.y_uc_dist_uc)./2;
sim.y_po_dist_do = (sim.sim.y_cu_dist_uc+sim.sim.y_uc_dist_cu)./2;

%% Fig. 6A-B. Model's Behavioral Bias
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
tl1.Layout.Tile = 3;
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
sim.y_val_po_dist_ds = (sim.y_val_cu_dist_cu);%+sim.y_val_uc_dist_uc);
sim.y_val_po_dist_do = (sim.y_val_cu_dist_uc);%+sim.y_val_uc_dist_cu);

x_axis = sim.x; 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
m_ds = mean(squeeze(sim.y_val_po_dist_ds(:,t_idx,:)),1);
m_do = mean(squeeze(sim.y_val_po_dist_do(:,t_idx,:)),1);

y_dog_ds = fit_curves(x_axis(x_idx), m_ds(x_idx), x_idx);
y_dog_do = fit_curves(x_axis(x_idx), m_do(x_idx), x_idx);


% figure
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 13;
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


