function [surf_lap] = eeg_laplaceTransform(EEG_ica)

addpath(genpath('/media/hdd/Sanchit/Exogenous_Project/Toolboxes/CSDtoolbox'));

% create Montage (M) 
for el=1:length(EEG_ica.chanlocs)
M.lab{el,1} = EEG_ica.chanlocs(el).labels;
end
M.theta = [EEG_ica.chanlocs(:).sph_theta]';
M.phi = [EEG_ica.chanlocs(:).sph_phi]';
locX = [EEG_ica.chanlocs(:).X]'; 
locY = [EEG_ica.chanlocs(:).Y]'; 
M.xy = [locX, locY];

[G,H] = GetGH(M);

surf_lap = CSD(EEG_ica.data, G, H);    


rmpath(genpath('/media/hdd/Sanchit/Exogenous_Project/Toolboxes/CSDtoolbox'));

end