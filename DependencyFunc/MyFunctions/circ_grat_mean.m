function m_deg = circ_grat_mean(x)
%%
addpath('/media/hdd/Sanchit/Exogenous_Project/Toolboxes/circstat-matlab-master'); 

% x is a vector/matrix of angles in degrees with range -90 to 90
% dim is the dimension across which the mean has to be taken
% output is in deg in the range -90 to 90
m_deg = nan(1,size(x,2)); 

for i = 1:size(x,2)
    
    idx = ~isnan(x(:,i)); 
    m_deg(1,i) = rad2deg( circ_mean(deg2rad(x(idx,i))*2, [], 1)/2 );
    
end

end