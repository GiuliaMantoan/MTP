% Procedure from Bai&Ng, JBES(2005), translated from GAUSS
% u is not demeaned!
function [hac,k] = nwX(u,prewhite,k)

% k is number of lags to be used
% if k is not provided or equals -1, it is determined automatically (Andrews 1991)
% if prewhite == 0 and k == 0, vcv-matrix equals standard vcv-matrix known from OLS 

n = rows(u);             % Number of observations
nreg = cols(u);          % Number of variables (moments/residuals)
rho = zeros(nreg,1);     % For storing AR(1) coefficients
sigma = zeros(nreg,1);   % For storing innovation variances
d = zeros(nreg,nreg);    % Identity used later in adjustment
beta = zeros(nreg,nreg); % To store VAR(1) coefficients if prewhitening

% Fit a VAR(1) model to each variable in u
if prewhite == 1
    v = zeros(n-1,nreg);
    reg = u(1:n-1,:);
    i = 1;
    while i<=nreg
        beta(:,i) = lscov( reg , u(2:n,i) );
        v(:,i) = u(2:n,i) - reg * beta(:,i);
        i = i + 1;
    end
else
    v = u;
end

if (nargin < 3)||(k==-1)
    i = 1;
    while i<=nreg
        rho(i) = lscov( v(1:rows(v)-1,i) , v(2:rows(v),i) );
        r = v(2:rows(v),i) - v(1:rows(v)-1,i) * rho(i);
        sigma(i) = r'*r / rows(v);
        i = i + 1;
    end


    bot = 0; top = 0;
    i = 1;
    while i<=nreg
        top = top + 4*(rho(i)^2)*sigma(i)^2 / (((1-rho(i))^6)*(1+rho(i))^2);
        bot = bot + sigma(i)^2 / ((1-rho(i))^4);
        i = i+1;
    end
    alpha = top/bot;
    k = ceil(1.1447*(alpha*n)^(1/3));
end

if k > n/2
    k = n/2;
end

%disp('truncation lag'); 
%disp(k);

% Initial term: standard covariance
vcv = v'*v / (n-1);
i = 1;
% Loop through lags and accumulate weighted autocovariances
while i<=k
    % x = i/k;
    w = 1 - i/(k+1);
    cov = v(i+1:rows(v),:)' * v(1:rows(v)-i,:) / (n-1);
    vcv = vcv + w*cov;
    cov = v(1:rows(v)-i,:)' * v(i+1:rows(v),:) / (n-1);
    vcv = vcv + w*cov;
    i = i + 1;
end
d = inv(eye(nreg)-beta'); % Transformation matrix
hac = d*vcv*d';

