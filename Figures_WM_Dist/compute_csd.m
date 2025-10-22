function [CSD] = compute_csd(X_resp,T_act)


X = wrap((X_resp*pi)/180);
T = wrap((T_act*pi)/180);

CSD = cstd(wrap(X-T)); 
return