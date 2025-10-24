function [x_axis, y_axis, x_val, y_val, area_bias] = WMA_biasAnalysis_bias(angleDiff, probeOrientation, distractorOrientation)

binWidth = 30;
x_val = -90:90;
y_val = nan(1,size(x_val,2)); 

y_axis = angleDiff;  

x_axis = circ_dist_grat(distractorOrientation, probeOrientation); 


iter = 0; 
for angle = x_val
    iter = iter+1;
    idx = abs(circ_dist_grat(x_axis, angle)) <=  binWidth/2; 
    y_val(iter) = nanmedian(y_axis(idx));
end
y_val = y_val - nanmean(y_val);

% Area under the bias curve 
idx = abs(x_val) <= 60 & abs(x_val) > 10; 
area_bias =  nanmean(y_val(idx).*sign(x_val(idx)));
end