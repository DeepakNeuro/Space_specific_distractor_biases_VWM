%{
modified pdist2() Matlab function with w_mat as an additional output
%}

function [D,W] = pdist2_mod(X, Y, varargin)
%PDIST2 Pairwise distance between two sets of observations.
%   D = PDIST2(X,Y) returns a matrix D containing the Euclidean distances
%   between each pair of observations in the MX-by-N data matrix X and
%   MY-by-N data matrix Y. Rows of X and Y correspond to observations,
%   and columns correspond to variables.
%
%   D = PDIST2(X,Y,DISTANCE) computes D using DISTANCE.  Choices are:
%
%       Mahalanobis distance, using the sample
%       covariance of X as computed by NANCOV.  To
%       compute the distance with a different
%       covariance, use
%       D = PDIST2(X,Y,'mahalanobis',C), where the
%       matrix C is symmetric and positive definite.


if nargin < 2
    error(message('stats:pdist2:TooFewInputs'));
end

[nx,p] = size(X);
[ny,py] = size(Y);
if py ~= p
    error(message('stats:pdist2:SizeMismatch'));
end

additionalArg = [];

% Custom covariance matrix
if ~isempty(varargin)
    arg = varargin{1};
    
    % Get the additional distance argument from the inputs
    if isnumeric(arg)
        additionalArg = arg;
    end
end

% Check data type
try
    outClass = superiorfloat(X,Y);
catch
    if isfloat(X)
        outClass = class(X);
    elseif isfloat(Y)
        outClass = class(Y);
    else
        outClass = 'double';
    end
end

if ~strcmp(class(X),outClass) || ~strcmp(class(Y),outClass)
    warning(message('stats:pdist2:DataConversion', outClass));
end
X = cast(X,outClass);
Y = cast(Y,outClass);
if  ~isreal(X) || ~isreal(Y)
    error(message('stats:pdist2:ComplexData'));
end

% Mahalanobis case in pdist2()
if isempty(additionalArg)
    if nx == 1
        error(message('stats:pdist2:tooFewXRowsForMah'));
    end
    additionalArg = nancov(X);
    [T,flag] = chol(additionalArg);
else %provide the covariance for mahalanobis
    if ~isequal(size(additionalArg),[p,p])
        error(message('stats:pdist2:InvalidCov'));
    end
    %cholcov will check whether the covariance is symmetric
    [T,flag] = cholcov(additionalArg,0);
end

if flag ~= 0
    error(message('stats:pdist2:InvalidCov'));
end

% Output variable
D = zeros(nx, ny, outClass);
W = zeros(nx, p, ny, outClass);

for i = 1:ny
    del = bsxfun(@minus, X, Y(i, :));
    dvc = del/T; % Displacement vectors
    dsq = sum(dvc.^2, 2);
    dsq = sqrt(dsq); % Distance array
    D(:, i) = dsq;
%     for k = 1:p
%         bv = zeros(size(dvc)); bv(:, k) = dvc(:, k); % Basis vector
%         pvm = sum(bv.*dvc, 2)./dsq; % Projected vector magnitutde
%         W(:, k, ny) = pvm; % pvm./dsq; % Weight
%     end
end