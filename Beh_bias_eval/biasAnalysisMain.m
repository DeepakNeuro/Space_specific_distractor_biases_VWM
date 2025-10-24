function [x_axis, y_axis, x_val, y_val, area_bias] = biasAnalysisMain(trialData, SubIdx)

% Calculate the raw angle difference
probeOrientation = [];
for i = 1:size(trialData,1)
    probeOrientation = [probeOrientation; trialData.initialAngles(i,trialData.WM_Probe(i))];     
end

distractorOrientation = trialData.distractorAngle; 
 
angleDiff = circ_dist_grat(trialData.ResponseWM, probeOrientation);

[x_axis, y_axis, x_val, y_val, area_bias] = WMA_biasAnalysis_bias(angleDiff, probeOrientation, distractorOrientation);
 
return;

