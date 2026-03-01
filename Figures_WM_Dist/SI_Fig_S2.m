a.tntu = load('../../Data/Behavior/t_nt_u_vonmisesfit_alldist_corrected.mat');
a.tu = load('../../Data/Behavior/t_u_vonmisesfit_alldist_corrected.mat');
a.tnt = load('../../Data/Behavior/t_nt_vonmisesfit_alldist.mat');
a.t = load('../../Data/Behavior/t_vonmisesfit_alldist.mat');


a.tntu.k = 4;%(2*2 + 3-1)
a.tu.k = 2;
a.t.k = 1;
a.tnt.k = 3;

for mix = ["tu", "tntu", "tnt", "t"]
    %% Compute AIC
    AIC_cu_ds.(sprintf('%s',mix)) = 2*a.(sprintf('%s',mix)).k - 2.*a.(sprintf('%s',mix)).LL_cu_ds;
    AIC_cu_do.(sprintf('%s',mix)) = 2*a.(sprintf('%s',mix)).k - 2.*a.(sprintf('%s',mix)).LL_cu_do;

    %% Compute AICc
    cfactor_ds = 2*((a.(sprintf('%s',mix)).k)^2 + a.(sprintf('%s',mix)).k)./(a.(sprintf('%s',mix)).n_obs_cu_ds - a.(sprintf('%s',mix)).k-1);
    cfactor_do = 2*((a.(sprintf('%s',mix)).k)^2 + a.(sprintf('%s',mix)).k)./(a.(sprintf('%s',mix)).n_obs_cu_do - a.(sprintf('%s',mix)).k-1);
    AICc_cu_ds.(sprintf('%s',mix)) = AIC_cu_ds.(sprintf('%s',mix))+cfactor_ds;
    AICc_cu_do.(sprintf('%s',mix)) = AIC_cu_do.(sprintf('%s',mix))+cfactor_do;


    %% Compute BIC
    BIC_cu_ds.(sprintf('%s',mix)) = log(a.(sprintf('%s',mix)).n_obs_cu_ds).*a.(sprintf('%s',mix)).k - 2.*a.(sprintf('%s',mix)).LL_cu_ds;
    BIC_cu_do.(sprintf('%s',mix)) = log(a.(sprintf('%s',mix)).n_obs_cu_do).*a.(sprintf('%s',mix)).k - 2.*a.(sprintf('%s',mix)).LL_cu_do;


end

%% plot for comparisiosn

%% Combine all plots into one tiled figure
figure;
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

%% ---------------- AICc comparison ----------------
m1 = AICc_cu_ds.tntu + AICc_cu_do.tntu;
m2 = AICc_cu_ds.tu   + AICc_cu_do.tu;
m3 = AICc_cu_ds.tnt  + AICc_cu_do.tnt;
m4 = AICc_cu_ds.t    + AICc_cu_do.t;

m = [m1', m2', m3', m4'];

nexttile;
violin_plot(m(:,1)-m(:,2));
title('AICc: t+nt+u vs t+u');
ylim([-10,25]);
ylabel('\Delta AICc (DS+DO)');

doPermutationTest(m(:,[1,2]),'both','both');

%% ---------------- BIC comparison ----------------
m1 = BIC_cu_ds.tntu + BIC_cu_do.tntu;
m2 = BIC_cu_ds.tu   + BIC_cu_do.tu;
m3 = BIC_cu_ds.tnt  + BIC_cu_do.tnt;
m4 = BIC_cu_ds.t    + BIC_cu_do.t;

m = [m1', m2', m3', m4'];

nexttile;
violin_plot(m(:,1)-m(:,2));
title('BIC: t+nt+u vs t+u');
ylim([-10,25]);
ylabel('\Delta BIC (DS+DO)');

doPermutationTest(m(:,[1,2]),'both','both');

%% ---------------- Alpha vs Beta scatter ----------------
nexttile;

pnt_cu_ds = a.tntu.B_cu_ds(4,:);
pt_cu_ds  = a.tntu.B_cu_ds(3,:);

m(:,1) = pnt_cu_ds;
m(:,2) = pt_cu_ds;

xlim([-0.05 1]);
ylim([-0.05 1]);
hold on;

r = refline(1);
r.Color = [1 1 1]*0.7;
r.LineStyle = '--';
r.LineWidth = 1;

scatter(m(:,1), m(:,2),80,'LineWidth',0.75, ...
    'MarkerEdgeColor', [1 1 1], ...
    'MarkerFaceColor', [0 0 0]);

xlabel('\beta');
ylabel('\alpha');
title('Von Mises mixture probs');

ax = gca;
ax.XTick = 0:0.2:1;
ax.YTick = ax.XTick;
axis square