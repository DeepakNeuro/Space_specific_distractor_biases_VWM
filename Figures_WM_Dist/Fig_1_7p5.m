%{
Code for reproducing the Fig. 1B,C,D,F,G & SI Fig. 1
%}

Figures_main;

%% Fig. 1B Histograms
m_cu = [beh.wm.raw.ang{1,1,:}]; % No distractor trials
m_cu = m_cu - nanmean(m_cu); 
m_cu = m_cu(:);
m_cu = deg2rad(m_cu(~isnan(m_cu)))*2;

m_uc = [beh.wm.raw.ang{2,1,:}]; % No distractor trials
m_uc = m_uc - nanmean(m_uc);
m_uc = m_uc(:);
m_uc = deg2rad(m_uc(~isnan(m_uc)))*2;

vm_fun = @(x,kappa, wG) (wG/2/pi + (1-wG)*1/2/pi/besseli(0,kappa)*exp(kappa*cos(x)));
ang = linspace(-pi,pi,100);
params = dmfit_wow(m_cu);
y_pdf_cu = vm_fun(ang,params(1), params(2));
params = dmfit_wow(m_uc);
y_pdf_uc = vm_fun(ang,params(1), params(2));

%figures
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 1;
tl1.Layout.TileSpan = [1 1];
tp = nexttile(tl1, 1);

histogram(m_cu, linspace(-pi,pi,37), 'Normalization', 'pdf', 'DisplayStyle','stairs','EdgeColor', color_map_att(1,:), 'LineWidth', 0.5);
hold on; 
histogram(m_uc, linspace(-pi,pi,37), 'Normalization', 'pdf', 'DisplayStyle','stairs','EdgeColor', color_map_att(2,:), 'LineWidth', 0.5);

plot(ang,y_pdf_cu,'LineWidth',2, 'Color', color_map_att(1,:)); 
plot(ang,y_pdf_uc,'LineWidth',2, 'Color', color_map_att(2,:)); 

ax = gca; 
ax.Box = 'off';
xlim([-1, 1]*pi);
ylim([0 0.8]); 
ax.XTick = linspace(-pi,pi,5);
ax.YTick = (0:0.4:0.8);
xlabel('Report Error (convert units)', 'FontSize', 8); 
ylabel('Prob. Desnity', 'FontSize', 8);

%% Fig. 1B Cued vs Uncued CSD (inset)
m = squeeze(beh.wm.f.JV10_error.csd(:,3,:))'; % No distractor trials (cued Vs uncued comparision)
pVal = doPermutationTest(m, 'smaller', 'left'); % left sided test

% figure
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 2;
tl1.Layout.TileSpan = [1 1];

tp = nexttile(tl1, 1);

xlim([10 60]);
ylim([10 60]);
r = refline(1); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
hold on; 
scatter(m(:,1), m(:,2),30,'LineWidth',0.75, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0 0 0]); 
xlabel('Cued'); 
ylabel('Uncued'); 
title('CSD (°)'); 
text(35,15,(sprintf('p = %0.3f', pVal))); 

ax = gca; 
ax.XTick = [10,60];
ax.YTick = ax.XTick; 
axis square

%% Fig. 1C CSD Violin Plots
m1 = squeeze(mean(beh.wm.f.JV10_error.csd, 1))'; % cued and uncued combined
m = [m1(:,3), m1(:,[1:2])];

% figure
tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 7;
tl1.Layout.TileSpan = [2 3];

tp = nexttile(tl1, 1);
xlim([0.65 3.35])
ylim([0, 70])
[vp, bx] = violin_plot(m, [color_map(6,:); color_map(1:2,:)]); 
ax = gca; 
ax.Box = 'off';
ax.XAxis.Visible = 'off';
ylim([0, 70])
ax.YTick = [0 35 70];
ylabel('CSD (°)')
tx = text([1,2,3], [0.1, 0.1, 0.1], {'No Dist.', 'Same', 'Opp.'}, 'HorizontalAlignment', 'center');
tx(1).Color = color_map(6,:); tx(2).Color = color_map(1,:); tx(3).Color = color_map(2,:); 

% ANOVA
d =  beh.wm.f.JV10_error.csd; %
clear dd; 
dd = d(:,[1 2 3], :);% cued and uncued distractor conditions (all distractor present trials)

g_att = repmat([1 1 1; 0 0 0],1,1,nSubjects);
g_dis = repmat([1 0 -1; 1 0 -1],1,1,nSubjects);
temp = permute([1:nSubjects]', [3,2,1]); 
g_sub = repmat(ones(2,3),1,1,nSubjects).*temp;

[p,tbl,stats] = anovan(dd(:),{g_att(:), g_dis(:), g_sub(:)},...
    'model','interaction', 'random', 3,...
    'varnames',{'Attention','DistLoc', 'Subject'});

% post-hoc signed-rank and BFs
doPermutationTest(m(:, [1,2]),'both','both'); % does signed-rank, t-test, BFs, permutaion tests
doPermutationTest(m(:, [1,3]),'both','both'); 
doPermutationTest(m(:, [2,3]),'both','both'); 


%% Fig. 1F,G Behavioral Bias plots (cued)
x_axis = beh.wm.bias.x.cu.ds(1,:); 
x_idx = (mod(x_axis, 15)==0) & (x_axis~=0);
x_axis_val = x_axis; x_axis_val(abs(x_axis_val)<10) = nan; 

% behavior bias cued
beh.wm.bias.y.cu.ds = beh.wm.bias.y.cu.ds;
beh.wm.bias.y.cu.do = beh.wm.bias.y.cu.do;

% mean bias
m_ds = mean(beh.wm.bias.y.cu.ds,1);
m_do = mean(beh.wm.bias.y.cu.do,1);

se_ds = (std(beh.wm.bias.y.cu.ds,[],1))./sqrt(nSubjects);
se_do = (std(beh.wm.bias.y.cu.do,[],1))./sqrt(nSubjects);

y_dog_ds = fit_curves(x_axis(x_idx), m_ds(x_idx), x_idx);
y_dog_do = fit_curves(x_axis(x_idx), m_do(x_idx), x_idx);

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 3;
tl1.Layout.TileSpan = [2, 2];
tp = nexttile(tl1,1);

hold on; 
ylim([-1, 1]*5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -10 20 20], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
p = errorbar(x_axis(x_idx), m_ds(x_idx), se_ds(x_idx), 'MarkerFaceColor', [1 1 1], 'Color', color_map(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0);
p.MarkerSize = 4; 
plot(x_axis, y_dog_ds, 'LineWidth', 2, 'Color', color_map(1,:)); 
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 5;
tl1.Layout.TileSpan = [2, 2];
tp = nexttile(tl1,1);

hold on; 
ylim([-1, 1]*5)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -10 20 20], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
p = errorbar(x_axis(x_idx), m_do(x_idx), se_do(x_idx), 'MarkerFaceColor', [1 1 1], 'Color', color_map(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0);
p.MarkerSize = 4; 
plot(x_axis, y_dog_do, 'LineWidth', 2, 'Color', color_map(2,:)); 
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-1:1:1]*3;


xlabel(tl1, 'Distractor - Stimulus Ori. (°)', 'FontSize', 8); 
ylabel(tl1, 'Response - Stimulus Ori. (°)', 'FontSize', 8); 
title(tl1, 'Orientation Bias', 'FontSize', 8)

%% 1F,G Behavioral Bias Violin Plots

% compute the area under bias for pooled (cued + uncued)
xx_cu = beh.wm.bias.x.cu.ds(1,:); 
yy_cu_ds = beh.wm.bias.y.cu.ds;
yy_cu_do = beh.wm.bias.y.cu.do;
idx = abs(xx_cu) <= 90 & abs(xx_cu) > 10;
beh.wm.bias.y.area.cu.ds =  nanmean(yy_cu_ds(:,idx).*sign(xx_cu(idx)),2);
beh.wm.bias.y.area.cu.do =  nanmean(yy_cu_do(:,idx).*sign(xx_cu(idx)),2);

% pooled
m_ds = beh.wm.bias.y.area.cu.ds; 
doPermutationTest(m_ds.*[1 0]); 
m_do = beh.wm.bias.y.area.cu.do; 
doPermutationTest(m_do.*[0 1]); 

tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 16;
tl1.Layout.TileSpan = [2 3];

tp = nexttile(tl1, 1);

ylim([-1, 1]*10)
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 0.75;
hold on; 
[vp, bp] = violin_plot(m_ds, color_map(1,:)); vp.ShowData = false; bp.LineWidth = 0.75;
ax = gca; 
ax.Box = 'off';
ax.XAxis.Visible = 'off';
ax.YTick = [-1:1:1]*4;
ylabel('Bias (°)')

tp = nexttile(tl1, 2);
ylim([-1, 1]*10)
hold on; 
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 0.75;
[vp, bp] = violin_plot(m_do, color_map(2,:)); vp.ShowData = false; bp.LineWidth = 0.75;
ax = gca; 
ax.Box = 'off';
ax.XAxis.Visible = 'off';
ax.YTick = [-1:1:1]*4;
ylabel('Bias (°)')



%% Fig 1 D (delta kappa DS - DO) comparisions
%script for fitting the von Mises to the pooled responses
Search_MLE_vonMisesFit_kappa_pooledErrors


% Plotting code
% Combine into arrays for plotting
means = [delta_k_ds_do_cu_mu, delta_k_ds_do_uc_mu];
errors = [delta_k_ds_do_cu_std, delta_k_ds_do_uc_std];

% Labels for each bar
labels = {'DS-DO_CU', 'DS-DO_UC'};

% Create horizontal bar plot
figure; hold on;
b = barh(means, 'FaceColor', [0.6 0.6 0.9], 'EdgeColor', 'none');

% Add horizontal error bars
for i = 1:numel(means)
    errorbar(means(i), i, errors(i), 'horizontal', 'k', 'LineStyle', 'none', 'LineWidth', 1.2);
end

% Customize axes
set(gca, 'YTick', 1:2, 'YTickLabel', labels, 'YDir', 'reverse'); % reverse for top-down order
xlabel('Δk value');
ylabel('Condition');
title('Horizontal Bar Plot with Error Bars');
box off;
grid on;
xlim padded;


%% SI figure S1
plotDefaults; 
tl = tiledlayout(8,6,'TileSpacing','compact', 'Padding', 'compact');
% uncued hemifield bias
xx_cu = beh.wm.bias.x.cu.ds(1,:); 
idx = abs(xx_cu) <= 90 & abs(xx_cu) > 10;
m_ds = mean(beh.wm.bias.y.uc.ds,1);
m_do = mean(beh.wm.bias.y.uc.do,1);

se_ds = (std(beh.wm.bias.y.uc.ds,[],1))./sqrt(nSubjects);
se_do = (std(beh.wm.bias.y.uc.do,[],1))./sqrt(nSubjects);

y_dog_ds = fit_curves(x_axis(x_idx), m_ds(x_idx), x_idx);
y_dog_do = fit_curves(x_axis(x_idx), m_do(x_idx), x_idx);

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 13;
tl1.Layout.TileSpan = [2, 2];
tp = nexttile(tl1,1);

hold on; 
ylim([-1, 1]*10)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -10 20 20], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
p = errorbar(x_axis(x_idx), m_ds(x_idx), se_ds(x_idx), 'MarkerFaceColor', [1 1 1], 'Color', color_map(1,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0);
p.MarkerSize = 4; 
plot(x_axis, y_dog_ds, 'LineWidth', 2, 'Color', color_map(1,:)); 
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-4:1:4]*2;

tl1 = tiledlayout(tl, 1, 1, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 15;
tl1.Layout.TileSpan = [2, 2];
tp = nexttile(tl1,1);

hold on; 
ylim([-1, 1]*6)
xlim([-1, 1]*90)
rect = rectangle('Position', [-10 -10 20 20], 'FaceColor', [0.9 0.9 0.9], 'LineStyle', 'none');
p = errorbar(x_axis(x_idx), m_do(x_idx), se_do(x_idx), 'MarkerFaceColor', [1 1 1], 'Color', color_map(2,:), 'LineStyle', 'none', 'Marker', 'o', 'LineWidth', 0.5, 'CapSize', 0);
p.MarkerSize = 4; 
plot(x_axis, y_dog_do, 'LineWidth', 2, 'Color', color_map(2,:)); 
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
r = vline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 1;
ax = gca; 
ax.Box = 'off';
ax.XTick = [-90:45:90];
ax.YTick = [-3:1:3]*2;

xlabel(tl1, 'Distractor - Stimulus Ori. (°)', 'FontSize', 8); 
ylabel(tl1, 'Response - Stimulus Ori. (°)', 'FontSize', 8); 
title(tl1, 'Orientation Bias', 'FontSize', 8)

%% violin plots SI fig S1
% uncued
m_ds = beh.wm.bias.y.area.uc.ds; 
doPermutationTest(m_ds.*[1 0],'larger','right'); 
m_do = beh.wm.bias.y.area.uc.do; 
doPermutationTest(m_do.*[1 0], 'smaller', 'left'); 

tl1 = tiledlayout(tl, 1, 2, 'TileSpacing','compact', 'Padding', 'compact');
tl1.Layout.Tile = 37;
tl1.Layout.TileSpan = [2 3];

tp = nexttile(tl1, 1);

ylim([-1, 1]*20)
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 0.75;
hold on; 
[vp, bp] = violin_plot(m_ds, color_map(1,:)); vp.ShowData = false; bp.LineWidth = 0.75;
ax = gca; 
ax.Box = 'off';
ax.XAxis.Visible = 'off';
ax.YTick = [-2:1:2]*4;
ylabel('Bias (°)')

tp = nexttile(tl1, 2);
ylim([-1, 1]*20)
hold on; 
r = hline(0); r.Color = [1 1 1]*0.7; r.LineStyle = '--'; r.LineWidth = 0.75;
[vp, bp] = violin_plot(m_do, color_map(2,:)); vp.ShowData = false; bp.LineWidth = 0.75;
ax = gca; 
ax.Box = 'off';
ax.XAxis.Visible = 'off';
ax.YTick = [-2:1:2]*4;
ylabel('Bias (°)')