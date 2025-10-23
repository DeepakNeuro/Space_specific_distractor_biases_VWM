function fitpars = dmfit_wow(x)

% Function by Priyanka
% x should be a vector from range -pi to pi

rng('default');

% find initial parameter estimate
maxLLH=-Inf;
kappa_vec = linspace(1,20,10);
w_vec = linspace(0+eps,1-eps,10);
    for jj = randperm(length(kappa_vec))
        kappa_hat = kappa_vec(jj);
        w_hat = w_vec(jj); 
        LLH = -mLLHfun([kappa_hat, w_hat],x);
        if LLH>maxLLH
            maxLLH=LLH;
            initpars = [kappa_hat, w_hat];
        else
            initpars = [kappa_hat, w_hat];
        end
    end
fitpars = fminsearch(@(pars) mLLHfun(pars,x), initpars);
return

function mLLH = mLLHfun(pars,x)
    wG = pars(2); % Fit uniform also
    kappa = pars(1);
    if kappa<0 || kappa>700
        mLLH = Inf;
    else
        mLLH = -sum(log(wG/2/pi + (1-wG)*1/2/pi/besseli(0,kappa)*exp(kappa*cos(x))));
    end
return;
