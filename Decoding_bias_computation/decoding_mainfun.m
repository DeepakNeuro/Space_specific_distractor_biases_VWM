function [dists_m, cos_m] = decoding_mainfun(dat, item_orientation)

%% Reduce dimensionality and zscore the data
dat = decodingTrc_pca(dat); 

addpath('/media/hdd/Sanchit/Exogenous_Project/Toolboxes/Decoding/Gokul_March_2021');
rng('default'); 

% setting up the orientation spaces to group the item orientations into
bin_width=180/6;
ang_steps=8; % 8 different orientation spaces to loop over
step_size = 180/16; 

angspace_temp= (-90:step_size:90)';
angspace_temp(end)=[];
for as=1:ang_steps
    angspaces(:,as)=angspace_temp+(as-1)*step_size/(ang_steps);
end


%% preallocate decoding output matrices for each item and epoch
dists_mem=nan(size(angspaces,1),ang_steps,size(dat,3), size(dat,1)); cos_mem=nan(ang_steps, size(dat,3), size(dat,1));

for a=1:ang_steps % loop through each orientation space
    
    %% decoding orientation
    [distance_cos,distances] = cog_mahal_compute(dat,item_orientation,angspaces(:,a), bin_width);
    
    dists_mem(:,a,:,:)= -permute(distances, [2 3 1]);
    cos_mem(a,:,:)= permute(distance_cos, [2 1]);
        
end
dists_m = dists_mem; %% preserving the angspaces dimension (Rm for normal)
cos_m = cos_mem; %% preserving the angspaces dimension (Rm for normal)
% dists_m = squeeze(mean(dists_mem, 2)); 
% cos_m = squeeze(mean(cos_mem, 1));
end