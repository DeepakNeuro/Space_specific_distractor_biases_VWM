function [vp, bp] = violin_plot(data, color_data)
% data should be nSamples x nGroups vector
% color_data should be a nGroups x 3 matrix

nSamples = size(data,1); 
nGroups = size(data,2); 
xGroup = zeros(nSamples,1) + (1:nGroups); 


if nargin < 2
    color_data = lines(nGroups); 
end

% Box Plot
bp = boxchart(xGroup(:), data(:), 'MarkerStyle', 'none', 'BoxFaceColor', [0 0 0], 'BoxFaceAlpha', 0, 'LineWidth', 1, 'BoxWidth', 0.1, ...
    'WhiskerLineColor', [0 0 0]); 


vp = violinplot(data); 

for i = 1:nGroups
    vp(i).ViolinColor = {color_data(i,:)}; vp(i).EdgeColor = color_data(i,:); 
    vp(i).ViolinAlpha = {[0.3]}; 
    vp(i).ShowData = false;
    vp(i).ScatterPlot.SizeData = 20; vp(i).ScatterPlot.MarkerEdgeColor = [1 1 1];
    vp(i).ShowBox = false; vp(i).ShowWhiskers = false; vp(i).ShowMedian = false;
end


end