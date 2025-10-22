function [fit_best] = fit_curves(x,y, x_idx, typeFit)

if nargin <=3
    typeFit = false; 
end

rng('default');

y_fit(1,:) = fit_dogAndUniform(x,  y);
% y_fit(2,:) = -(fit_dogAndUniform(x, -y));
% [params, y_fit(3,:)] = fit_sin(x,  y);
% 
% if params(2) < 75 && ~typeFit
%     y_fit(3,:) = 0; % Don't fit a sin that has period less than 75 deg
% end

% if typeFit
%     y_fit(typeFit,:) = 0; % Don't fit a particular type of curve -- replace with zeros
% end

% r_fits = corr(y', y_fit(:, x_idx)');
% [~, max_idx] = max(r_fits); 

fit_best = y_fit(1,:); 
% fit_best = y_fit(3,:); 


end
