function r_deg = circ_dist_grat(x,y)

% x and y are angles in degrees with range -90 to 90
% r_deg is also in degree with the range -90 to 90

r = angle(exp(1i*deg2rad(2*x))./exp(1i*deg2rad(2*y)));

r_deg = rad2deg(r)/2; 

end