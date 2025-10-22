function [fig, p1, r] = plt_scatter_xy(dd)

% fig = figure('Position', [1 1 300 300]); % nan; %
p1 = scatter(dd(:,1), dd(:,2), 'filled');

ax = gca;
% ax.LineWidth = 1.5;
% ax.FontSize = 9;

ax.Box = 'off';

% axis tight; 
limMin = min(ax.XLim(1), ax.YLim(1));
limMax = max(ax.XLim(2), ax.YLim(2));

ax.XLim = [limMin, limMax];
ax.YLim = [limMin, limMax];
if length(ax.XTick) <= length(ax.YTick)
    ax.YTick = ax.XTick;
else
    ax.XTick = ax.YTick;
end

r = refline(1);
r.Color = [1 1 1]*0.7;
r.LineStyle = '--';
r.LineWidth = 0.7;

% axis equal;
axis square;

end