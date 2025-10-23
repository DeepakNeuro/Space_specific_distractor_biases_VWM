function [B, LL, W] = mixtureFunction (X, T, NT, B0)
% mixtureFunction (X, T, NT, B0) 
%   Expectation Maximization function called by mixtureFit with starting parameters specified in B0.
%
%   Paul Bays | bayslab.com | Licence GPL-2.0 | 2017-02-20

if (nargin<2 || size(X,2)>1 || size(T,2)>1 || size(X,1)~=size(T,1) || nargin>2 && ~isempty(NT) && (size(NT,1)~=size(X,1) || size(NT,1)~=size(T,1)))     
    error('Input is not correctly dimensioned');
    return; 
end

if (nargin>3 && (B0(1)<0 || B0(2)<0 || any(B0(3:5)<0) || any(B0(3:5)>1) || abs(sum(B0(3:5))-1) > 10^-6))
    error('Invalid model parameters');
    return;
end

MaxIter = 10^4; MaxdLL = 10^-4;

n = size(X,1); 

if (nargin<3) 
    NT = zeros(n,0); nn = 0;
else
    nn = size(NT,2);
end

% Default starting parameters
if (nargin<4)    
    K1 = 5; K2 = 5; Pt = 0.5; 
    if (nn>0) Pn = 0.3; else Pn = 0; end
    Pu = 1-Pt-Pn;
else
    K1 = B0(1);
    K2 = B0(2);
    Pt = B0(3); Pn = B0(4); Pu = B0(5);
end


%% previous lines (accepts already convered to -pi to pi and perfoms wrapping around)
% E  = X-T; E = mod(E + pi, 2*pi) - pi;
% NE = repmat(X,1,nn)-NT;  NE = mod(NE + pi, 2*pi) - pi;

%% modified by Deepak Raya for -90 to 90 deg to -pi to pi rad (main thing keep it check thsese carefully)
% E  = deg2rad(X-T)*2; E = E(~isnan(X)); E = mod(E + pi, 2*pi) - pi;
% NE = deg2rad(repmat(X,1,nn)-NT)*2;  NE = NE(~isnan(X)); NE = mod(NE + pi, 2*pi) - pi;

%% modified by Deepak Raya for permtute kappa analysis
E  = X-T; E = mod(E + pi, 2*pi) - pi;
NE = repmat(X,1,nn)-NT; NE = mod(NE + pi, 2*pi) - pi;



n = size(E,1);

LL = nan; dLL = nan; iter = 0;

while (1)
    iter = iter + 1;
    
    Wt = Pt * vonmisespdf(E,0,K1);
    Wu = Pu * ones(n,1)/(2*pi);

    if nn==0
        Wn = zeros(size(NE));
    else
        Wn = Pn/nn * vonmisespdf(NE,0,K2);
    end
    
    W = sum([Wt Wn Wu],2);
    
    dLL = LL-sum(log(W));
    LL = sum(log(W));
    if (abs(dLL) < MaxdLL | iter > MaxIter) break; end
    
    Pt = sum(Wt./W)/n;
    Pn = sum(sum(Wn,2)./W)/n; 
    Pu = sum(Wu./W)/n;
            
    rw1 = (Wt./W);
    rw2 = (Wn./repmat(W,1,nn));
    
    S1 = sin(E); C1 = cos(E);
    S2 = sin(NE); C2 = cos(NE);
    
    r1 =[sum(sum(S1.*rw1)) sum(sum(C1.*rw1))];
    r2 =[sum(sum(S2.*rw2)) sum(sum(C2.*rw2))];
    
    if sum(sum(rw1))==0
        K1 = 0;
    else
        R1 = sqrt(sum(r1.^2))/sum(sum(rw1));        
        K1 = fastA1inv(R1);
    end
    
    if sum(sum(rw2))==0
        K2 = 0;
    else
        R2 = sqrt(sum(r2.^2))/sum(sum(rw2));        
        K2 = fastA1inv(R2);
    end  
    
    if n<=15
        if K1<2
            K1 = max(K1-2/(n*K1), 0);
        else
            K1 = K1 * (n-1)^3/(n^3+n);
        end
    end
    
    if n<=15
        if K2<2
            K2 = max(K2-2/(n*K2), 0);
        else
            K2 = K2 * (n-1)^3/(n^3+n);
        end
    end  
end

if iter>MaxIter
    warning('mixtureFunction:MaxIter','Maximum iteration limit exceeded.');
    B = [NaN NaN NaN NaN NaN]; LL = NaN; W = NaN;
else  
    B = [K1 K2 Pt Pn Pu]; W = [Wt sum(Wn,2) Wu]./sum([Wt Wn Wu],2);
end