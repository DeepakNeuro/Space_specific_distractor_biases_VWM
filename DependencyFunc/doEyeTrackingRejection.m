function [trialDataRejected, ETrejtIdx] =  doEyeTrackingRejection(subjectNumber, trialData)

BadETdataSubjects = [113, 115, 121, 126, 128, 133]; % Data was not saved for #113; All other had high rejection rates

trialDataRejected = trialData;
ETrejtIdx = false(size(trialData.ResponseWM));

if ~ismember(subjectNumber, BadETdataSubjects)
    %% Reject based on Eye movements
    fileET_WM = dir(sprintf('/media/hdd/Sanchit/Exogenous_Project/Analysis Data/ET/Sub_%d_ET_reject_20210720.mat',subjectNumber));
    load(sprintf('%s/%s',fileET_WM.folder,fileET_WM.name),'ETreject');
    
    %% Choose trials that are to be rejected
    
    % Fixation reject
    [~, side] = min([length(ETreject.Overall.L), length(ETreject.Overall.R)]);
    
    if side==1
        reject_idx = ETreject.Overall.L;
    elseif side==2
        reject_idx = ETreject.Overall.R;
    end
    
    %% Reject trials
    trialDataRejected.ResponseWM(reject_idx) = nan;
    trialDataRejected.RT_WM(reject_idx,:) = nan;
    
    ETrejtIdx(reject_idx) = true;
else
    fprintf('\nNo ET rejection done for Sub%d', subjectNumber); 
end

end