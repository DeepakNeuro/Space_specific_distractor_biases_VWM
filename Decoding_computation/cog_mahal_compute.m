function  [cos_amp, d_tune, cos_tune, w_mat] = cog_mahal_compute(data, theta, angspace, bin_width)


%% Input
% data     : EEG data (trial x channels x samples)
% theta    : Trial angles
% angspace : The angular space of the tuning curve, where each number is the center of the angle bins relative to the orienation of the test-trial (1 x bins)
% bin_width: The range of orientations that are included in each angle bin

%% Output
% cos_amp : The cosine amplitude of the tuning curve. A summary value for the decoding accuracy (matrix - trials x samples)
% d_tune  : The trial-wise distance tuning curves (matrix - trials x bins x samples)
% weights : Weights of electrodes in decoding (matrix - electrodes x samples)


%% Authors and Credits
% Cosine projection - (see Sprague, T. C., Ester, E. F. & Serences, J. T. Restoring Latent Visual Working Memory Representations in Human Cortex. Neuron 91, 694–707 (2016))
% This code was written by Michael Wolff and Mark Stokes (see Wolff et al., 2017, Nature Neuroscience)
% modified by Gokul T Anandan (gokul.t.anandan@gmail.com)
% modified by Sanchit Gupta (sanchitgupta@iisc.ac.in)
% modified by Deepak Raya (deepakvr@iisc.ac.in)

%% Execution
if size(angspace,2) == length(angspace)
    angspace = angspace';
end

% Convert degrees to radians
theta     = theta*2*pi/360;
angspace  = angspace*2*pi/360;
bin_width = bin_width*2*pi/360;

% Initialize variables
d_tune   = nan(size(data, 1), length(angspace), size(data, 3)); % trials x bins x samples
cos_amp  = nan(size(data, 1), size(data, 3)); % trials x samples
cos_tune = nan(size(d_tune)); % trials x bins x samples
trl_ind  = 1:size(data, 1); % trial indices
w_mat    = zeros(length(angspace), size(data, 2), size(data, 1), size(data, 3)); % bins x electrodes x trials x samples

reverseStr = ''; % This is simply to present the percentage completed

for trl = 1:size(data,1) % Trial number
    trn_dat = data(setdiff(trl_ind ,trl), :, :); % Training data
    trn_angle = theta(setdiff(trl_ind, trl)); % Traning angles
    
    m = nan(length(angspace), size(trn_dat,2), size(trn_dat,3));
    for b = 1:length(angspace)
        % Average the training data into orientation bins relative to the test-trial's orientation
        m(b, :, :) = mean(trn_dat(abs(angle(exp(1i*(2*trn_angle))./exp(1i*(2*(theta(trl)-angspace(b))))))/2 < bin_width, :, :), 1);
    end
    
    msg = sprintf('%d percent', round((trl/size(data, 1))*100));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    for ti = 1:size(data, 3)
        if ~isnan(trn_dat(:, :, ti))
            % The covariance matrix is computed for each time-point and excluding the test-trial
            sigma = covdiag(trn_dat(:, :, ti)); % here covariance is estimated using shrinkage (see below from Ledoit & Wolf, 2004).
            
            % Calculates the distances between the trial and all angle bins
            [d_tune(trl, :, ti), ~] = pdist2_mod(squeeze(m(:, :, ti)), squeeze(data(trl, :, ti)), sigma); % modified pdist2() with w_mat as an additional output
        end
    end
end

% De-mean the tuning curves along the bin dimension
d_tune = d_tune - repmat(nanmean(d_tune, 2), [1 size(d_tune, 2) 1]);

% Evaluate cosine average of angspace for tuning curves from mean-centered d_tune
% Since the decoding tuning curve should resemble a reversed cosine (higher distance = higher value), the value is reversed for ease of interpretation, so that high = high decoding
for trl = 1:size(d_tune, 1)
    for ti = 1:size(d_tune, 3)
        cos_amp(trl,ti) = -(mean(cos(angspace).*squeeze(d_tune(trl,:,ti))'));
    end
end

end