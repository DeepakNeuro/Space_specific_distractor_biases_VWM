function [ang, y_pdf, k, w] = fit_vM_simple(m)


vm_fun = @(x,kappa, wG) (wG/2/pi + (1-wG)*1/2/pi/besseli(0,kappa)*exp(kappa*cos(x)));

ang = linspace(-pi,pi,100);
params = dmfit_wow(m);
y_pdf = vm_fun(ang,params(1), params(2));

k = params(1);
w = params(2);

return