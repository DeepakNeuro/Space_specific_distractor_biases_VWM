% dog  Fit first derivative of Gaussian
%
% [ALPHA,SIGMA, AMP] = DOG(X,Y) fits first derivative of Gaussian to
% x,y-data by minimizing the sum of squared residuals. The output parameter
% ALPHA controls amplitude and SIGMA is the standard deviation of the
% Gaussian distribution and controls width of the resulting curve, given by
% y = normpdf(x,0,SIGMA).*x.*ALPHA. AMP is the peak amplitude. 
% 
%   Example
%       x = -50:1:50;
%       a = 40;
%       s = 15;
%       Fn = @(x,s) (1./sqrt(2.*pi.*sigma.^2)).*exp(-(x.^2)./(2.*sigma.^2));
%       y = Fn(x,s).*x.*a+randn(1,length(x)).*3;
%       [alpha,sigma] = dog(x,y);
%       ydog = Fn(x,sigma).*x.*alpha;
%       plot(x,y,'Displayname','Data');hold on;
%       plot(x,ydog,'Displayname','DoG'); legend
%   
%   Modified example which will work -- the example given above may not work
%   because many variable names coincide with function names. 
%
%   x = -50:1:50;
%   a = 40;
%   s = 15;
%   Fn = @(x,s) (1./sqrt(2.*pi.*s.^2)).*exp(-(x.^2)./(2.*s.^2));
%   y = Fn(x,s).*x.*a+randn(1,length(x)).*3;
%   [al,s, amp] = dog(x,y);
%   ydog = Fn(x,s).*x.*al;
%   plot(x,y,'Displayname','Data');hold on;
%   plot(x,ydog,'Displayname','DoG'); legend
%
% Taken from: https://www.mathworks.com/matlabcentral/fileexchange/70203-dog-fit-first-derivative-of-gaussian?s_tid=prof_contriblnk 

function [ydog, mu, s, a, b, amp] = fit_dogAndUniform(x,y)
rng('default');

Fn = @(x,mu,s) (1./sqrt(2.*pi.*s.^2)).*exp(-((x - mu).^2)./(2.*s.^2));

F = @(pr) sum((Fn(x,pr(1),pr(2)).*(x-pr(1)).*pr(3)  + pr(4)- y).^2);
 
opt = fmincon(F,[0,30,1,0],[],[],[],[],[-10,0+eps,-120,-60],[10,90,120,60]);

[mu,s,a,b] = deal(opt(1),opt(2),opt(3),opt(4));
xi = linspace(-90,90,181);
amp = max(Fn(xi,mu,s).*(xi-mu).*a + b);

ydog = Fn(xi,mu,s).*(xi-mu).*a  + b;
end
