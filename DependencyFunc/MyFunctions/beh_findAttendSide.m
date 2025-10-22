function [idxAttendLeft idxAttendRight] = beh_findAttendSide(subjectNumber, trialData)


idxAttendRight = trialData.WM_cue==2;
idxAttendLeft = trialData.WM_cue==1;


end