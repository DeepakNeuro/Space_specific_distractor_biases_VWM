%{
Code for reproducing SI figure coressponding to the ominibus analyis. 
%}
Figures_main

%% Load Data
Behcorr.time = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Comb_DistractorProject/Comb_March_24_2023_09_00_55.mat');
Behcorr.decoding = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Comb_DistractorProject/Comb_May_21_2024_16_17_19.mat');
Behcorr.gating = load('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/Comb_DistractorProject/Comb_May_21_2024_15_27_19.mat');

%%
m_lo(1, :,:) = Behcorr.time.comb.lo.wm.bias.y.cu.ds(2:24, :);
m_hi(1, :,:) = Behcorr.time.comb.hi.wm.bias.y.cu.ds(2:24, :);

m_lo(2, :,:) = Behcorr.decoding.comb.lo.wm.bias.y.cu.ds;
m_hi(2, :,:) = Behcorr.decoding.comb.hi.wm.bias.y.cu.ds;

m_lo(3, :,:) = Behcorr.gating.comb.hi.wm.bias.y.cu.ds; 
m_hi(3, :,:) = Behcorr.gating.comb.lo.wm.bias.y.cu.ds; 

m_lo_pool = cat(1, squeeze(m_lo(1, :,:)),squeeze(m_lo(2, :,:)),squeeze(m_lo(3, :,:)));
m_hi_pool = cat(1, squeeze(m_hi(1, :,:)),squeeze(m_hi(2, :,:)),squeeze(m_hi(3, :,:)));


%%
m_pool = [Behcorr.time.comb.lo.wm.bias.y.area.cu.ds(2:24), Behcorr.time.comb.hi.wm.bias.y.area.cu.ds(2:24);
Behcorr.gating.comb.hi.wm.bias.y.area.cu.ds, Behcorr.gating.comb.lo.wm.bias.y.area.cu.ds;
Behcorr.decoding.comb.lo.wm.bias.y.area.cu.ds, Behcorr.decoding.comb.hi.wm.bias.y.area.cu.ds];


%% pooled plots
x_axis = Behcorr.decoding.comb.lo.wm.bias.x.cu.ds(1,:); 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);


m_lo_plot = nanmean(m_lo_pool,1);
m_hi_plot  = nanmean(m_hi_pool,1);

se_lo = (nanstd(m_lo_pool,[],1))./(sqrt(nEEG*3));
se_hi = (nanstd(m_hi_pool,[],1))./(sqrt(nEEG*3));

k_lo_plot  = warp_mov_mean(m_lo_plot(x_idx),3);
k_hi_plot  = warp_mov_mean(m_hi_plot(x_idx),3);

y_dog_lo_plot  = fit_curves(x_axis(x_idx), k_lo_plot , x_idx); 
y_dog_hi_plot  = fit_curves(x_axis(x_idx), k_hi_plot , x_idx);

plotDefaults; 
tl = tiledlayout(8,6,'TileSpacing','compact', 'Padding', 'compact');

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);
hold on; 
ylim([-1, 1]*6)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -20 20 40], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
p = errorbar(x_axis(x_idx), k_lo_plot, se_lo(x_idx), 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4); hold on;
p = errorbar(x_axis(x_idx), k_hi_plot, se_hi(x_idx), 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4); hold on;
plot(x_axis, y_dog_lo_plot, 'LineWidth', 2, 'Color', color_map_split(1,:)); hold on;
plot(x_axis, y_dog_hi_plot, 'LineWidth', 2, 'Color', color_map_split(2,:)); hold on;
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*4;
title('Time');
tx = text([15, 15],[-4, -5],{'Early/Strong Gat.','Late/Weak Gat.'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Actual Ori. (°)')
ylabel('Response - Actual Ori. (°)')

% %% permutation test
doPermutationTest(m_pool,'smaller','left');
doPermutationTest(m_pool(:,1)*[1 0],'larger','right');
doPermutationTest(m_pool(:,2)*[1 0],'larger','right');
%% violin plot
tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 4;
tl1.Layout.TileSpan = [2 3];
tp = nexttile(tl1, 1);

xlim([0.5 2.5])
ylim([-12, 14]);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
hold on;
[vp, bp] = violin_plot(m_pool, color_map_split); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [-1:1:1]*10;
ax.XAxis.Visible = 'off';
ylabel('Bias (°)'); 
% tx = text([1,2], [12, 12], {'n.s.', 'n.s.'}, 'HorizontalAlignment', 'center');
tx(2).FontSize = 12;
tx(1).FontSize = 12; 
% tx = text([1,2], [-12, -12], {'Early/Weak Enco./Strong Gat.','Late/Strong Enco./Weak Gat.'}, 'HorizontalAlignment', 'center');
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);




%%
m_lo(1, :,:) = Behcorr.time.comb.lo.wm.bias.y.cu.do(2:24, :);
m_hi(1, :,:) = Behcorr.time.comb.hi.wm.bias.y.cu.do(2:24, :);

m_lo(2, :,:) = Behcorr.decoding.comb.lo.wm.bias.y.cu.do;
m_hi(2, :,:) = Behcorr.decoding.comb.hi.wm.bias.y.cu.do;

m_lo(3, :,:) = Behcorr.gating.comb.hi.wm.bias.y.cu.do; 
m_hi(3, :,:) = Behcorr.gating.comb.lo.wm.bias.y.cu.do; 

m_lo_pool = cat(1, squeeze(m_lo(1, :,:)),squeeze(m_lo(2, :,:)),squeeze(m_lo(3, :,:)));
m_hi_pool = cat(1, squeeze(m_hi(1, :,:)),squeeze(m_hi(2, :,:)),squeeze(m_hi(3, :,:)));


%%
m_pool = [Behcorr.time.comb.lo.wm.bias.y.area.cu.do(2:24), Behcorr.time.comb.hi.wm.bias.y.area.cu.do(2:24);
Behcorr.gating.comb.hi.wm.bias.y.area.cu.do, Behcorr.gating.comb.lo.wm.bias.y.area.cu.do;
Behcorr.decoding.comb.lo.wm.bias.y.area.cu.do, Behcorr.decoding.comb.hi.wm.bias.y.area.cu.do];


%% pooled plots
x_axis = Behcorr.decoding.comb.lo.wm.bias.x.cu.ds(1,:); 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);


m_lo_plot = nanmean(m_lo_pool,1);
m_hi_plot  = nanmean(m_hi_pool,1);

se_lo = (nanstd(m_lo_pool,[],1))./(sqrt(nEEG*3));
se_hi = (nanstd(m_hi_pool,[],1))./(sqrt(nEEG*3));

k_lo_plot  = warp_mov_mean(m_lo_plot(x_idx),3);
k_hi_plot  = warp_mov_mean(m_hi_plot(x_idx),3);

y_dog_lo_plot  = fit_curves(x_axis(x_idx), k_lo_plot , x_idx); 
y_dog_hi_plot  = fit_curves(x_axis(x_idx), k_hi_plot , x_idx);

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 13;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);
hold on; 
ylim([-1, 1]*6)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -20 20 40], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
p = errorbar(x_axis(x_idx), k_lo_plot, se_lo(x_idx), 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4); hold on;
p = errorbar(x_axis(x_idx), k_hi_plot, se_hi(x_idx), 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4); hold on;
plot(x_axis, y_dog_lo_plot, 'LineWidth', 2, 'Color', color_map_split(1,:)); hold on;
plot(x_axis, y_dog_hi_plot, 'LineWidth', 2, 'Color', color_map_split(2,:)); hold on;
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*4;
title('Time');
tx = text([15, 15],[-4, -5],{'Early/Strong Gat.','Late/Weak Gat.'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Actual Ori. (°)')
ylabel('Response - Actual Ori. (°)')

% %% permutation test
doPermutationTest(m_pool,'larger','right');
doPermutationTest(m_pool(:,1)*[1 0],'smaller','left');
doPermutationTest(m_pool(:,2)*[1 0],'smaller','left');
%% violin plot
tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 15;
tl1.Layout.TileSpan = [2 3];
tp = nexttile(tl1, 1);

xlim([0.5 2.5])
ylim([-12, 14]);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
hold on;
[vp, bp] = violin_plot(m_pool, color_map_split); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [-1:1:1]*10;
ax.XAxis.Visible = 'off';
ylabel('Bias (°)'); 
% tx = text([1,2], [12, 12], {'n.s.', 'n.s.'}, 'HorizontalAlignment', 'center');
tx(2).FontSize = 12;
tx(1).FontSize = 12; 
% tx = text([1,2], [-12, -12], {'Early/Weak Enco./Strong Gat.','Late/Strong Enco./Weak Gat.'}, 'HorizontalAlignment', 'center');
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);


%%  %%%%%%
%%
m_lo(1, :,:) = Behcorr.time.comb.lo.wm.bias.y.cu.ds(2:24, :) - Behcorr.time.comb.lo.wm.bias.y.cu.do(2:24, :);
m_hi(1, :,:) = Behcorr.time.comb.hi.wm.bias.y.cu.ds(2:24, :) - Behcorr.time.comb.lo.wm.bias.y.cu.do(2:24, :);

% m_lo(2, :,:) = Behcorr.decoding.comb.lo.wm.bias.y.cu.ds - Behcorr.decoding.comb.lo.wm.bias.y.cu.do;
% m_hi(2, :,:) = Behcorr.decoding.comb.hi.wm.bias.y.cu.ds - Behcorr.decoding.comb.hi.wm.bias.y.cu.do;

m_lo(2, :,:) = Behcorr.gating.comb.hi.wm.bias.y.cu.ds - Behcorr.gating.comb.hi.wm.bias.y.cu.do; 
m_hi(2, :,:) = Behcorr.gating.comb.lo.wm.bias.y.cu.ds - Behcorr.gating.comb.lo.wm.bias.y.cu.do; 

m_lo_pool = cat(1, squeeze(m_lo(1, :,:)),squeeze(m_lo(2, :,:)));%,squeeze(m_lo(3, :,:)));
m_hi_pool = cat(1, squeeze(m_hi(1, :,:)),squeeze(m_hi(2, :,:)));%,squeeze(m_hi(3, :,:)));


%%
m_pool = [Behcorr.time.comb.lo.wm.bias.y.area.cu.ds(2:24) - Behcorr.time.comb.lo.wm.bias.y.area.cu.do(2:24),...
         Behcorr.time.comb.hi.wm.bias.y.area.cu.ds(2:24) - Behcorr.time.comb.hi.wm.bias.y.area.cu.do(2:24);...
         Behcorr.gating.comb.hi.wm.bias.y.area.cu.ds - Behcorr.gating.comb.hi.wm.bias.y.area.cu.do, ...
         Behcorr.gating.comb.lo.wm.bias.y.area.cu.ds -  Behcorr.gating.comb.lo.wm.bias.y.area.cu.do;...
         Behcorr.decoding.comb.lo.wm.bias.y.area.cu.ds - Behcorr.decoding.comb.lo.wm.bias.y.area.cu.do,...
         Behcorr.decoding.comb.hi.wm.bias.y.area.cu.ds - Behcorr.decoding.comb.hi.wm.bias.y.area.cu.do];

%% prep data from rmANOVA
data = m_pool(:);
group = [repmat({'Group1'},46,1); repmat({'Group2'},46,1)];

% Run one-way ANOVA
[p, tbl, stats] = anova1(data, group);

% Optional: post-hoc multiple comparisons
multcompare(stats);



%% pooled plots
x_axis = Behcorr.decoding.comb.lo.wm.bias.x.cu.ds(1,:); 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);


m_lo_plot = nanmean(m_lo_pool,1);
m_hi_plot  = nanmean(m_hi_pool,1);

se_lo = (nanstd(m_lo_pool,[],1))./(sqrt(nEEG*3));
se_hi = (nanstd(m_hi_pool,[],1))./(sqrt(nEEG*3));

k_lo_plot  = warp_mov_mean(m_lo_plot(x_idx),3);
k_hi_plot  = warp_mov_mean(m_hi_plot(x_idx),3);

y_dog_lo_plot  = fit_curves(x_axis(x_idx), k_lo_plot , x_idx); 
y_dog_hi_plot  = fit_curves(x_axis(x_idx), k_hi_plot , x_idx);

plotDefaults; 
tl = tiledlayout(8,6,'TileSpacing','compact', 'Padding', 'compact');

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [2, 2];

tp = nexttile(tl1, 1);
hold on; 
ylim([-1, 1]*6)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -20 20 40], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
p = errorbar(x_axis(x_idx), k_lo_plot, se_lo(x_idx), 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4); hold on;
p = errorbar(x_axis(x_idx), k_hi_plot, se_hi(x_idx), 'MarkerFaceColor', [1 1 1], 'Color', color_map_split(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0, 'MarkerSize', 4); hold on;
plot(x_axis, y_dog_lo_plot, 'LineWidth', 2, 'Color', color_map_split(1,:)); hold on;
plot(x_axis, y_dog_hi_plot, 'LineWidth', 2, 'Color', color_map_split(2,:)); hold on;
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*4;
title('Time');
tx = text([15, 15],[-4, -5],{'Early/Strong Gat.','Late/Weak Gat.'}); 
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
xlabel('Distractor - Actual Ori. (°)')
ylabel('Response - Actual Ori. (°)')

% %% permutation test
 
doPermutationTest(m_pool,'smaller','left');
doPermutationTest(m_pool(:,1)*[1 0],'larger','right');
doPermutationTest(m_pool(:,2)*[1 0],'larger','right');
%% violin plot
tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 4;
tl1.Layout.TileSpan = [2 3];
tp = nexttile(tl1, 1);

xlim([0.5 2.5])
ylim([-12, 14]);
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
hold on;
[vp, bp] = violin_plot(m_pool, color_map_split); 
ax = gca; 
ax.Box = 'off';
ax.YTick = [-1:1:1]*10;
ax.XAxis.Visible = 'off';
ylabel('Bias (°)'); 
% tx = text([1,2], [12, 12], {'n.s.', 'n.s.'}, 'HorizontalAlignment', 'center');
tx(2).FontSize = 12;
tx(1).FontSize = 12; 
% tx = text([1,2], [-12, -12], {'Early/Weak Enco./Strong Gat.','Late/Strong Enco./Weak Gat.'}, 'HorizontalAlignment', 'center');
tx(1).Color = color_map_split(1,:); tx(2).Color = color_map_split(2,:);
