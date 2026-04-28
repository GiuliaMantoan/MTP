function [fittedParams,fittedQuantiles] = fit_tpn_to_quantiles_SM(data,quantiles,mode0, sigma1_0, sigma2_0)

% INPUT:
% data: a 1 x M vector of fitted quantiles
% quantiles: the actual quantiles. For example: [0.05 ... 0.5 ... 0.95]
% mode0, sigma_10 and sigma_20: initial guesses for the parameters

% OUTPT
% fittedParams: parameters from the 2-piece normal 
% fittedQuantiles: quantiles fitted from the minimization pb 

data      = data(:)'; % check that data are in rows
data = sort(data); % sort the data

% select the quantiles of interest
[~,jq50] = min(abs(quantiles-.50)); % median
[~,jq25] = min(abs(quantiles-.25)); % lower qt
[~,jq75] = min(abs(quantiles-.75)); % upper qt

% Initial guess based on empirical quantiles
if nargin < 3

    mode0    = data(jq50); % mode
    
    z75      = norminv(0.75, 0, 1);
    z25      = norminv(0.25, 0, 1);

    sigma1_0 = (data(jq50) - data(jq25)) / abs(z25); % s1
    sigma2_0 = (data(jq75) - data(jq50)) / z75; %s2

end

paramsInit = [mode0, sigma1_0, sigma2_0]; % initial guess

% Bounds: ensure sigmas positive
LB = [-Inf, 1e-5 , 1e-5];
UB = [ Inf, Inf, Inf];

opts = optimset('Display','off','Algorithm','interior-point', ...
                   'MaxFunEvals',1e4,'MaxIter',1e4,'TolX',1e-6,'TolCon',1e-6); % option for minimzation


resid = @(x) ( two_part_normal_inverse_cdf(quantiles([jq25, jq50, jq75]), x(1), x(2), x(3)) - data([jq25, jq50, jq75]) ); % residual from actual qt and qt obtained with the parameters
objfun = @(x) norm(resid(x)); % L2 norm 

sol = fmincon(objfun, paramsInit, [], [], [], [], LB, UB, [], opts ); % get the parameters that give you the min residual

% temporary fix for those cases in which we obtain a value for the std
% deviations that is very small
% if that's the case, we fit a normal distribution 
% if both sigmas are small then thre are some problems and so worth
% checking the data at disposal
% UPDATED: in accordance with Giulia, if one of the two sd is small we set
% it to 0.01
if sol(2) < 1e-2 && sol(3) > 1e-2
    sol(2) = 0.01;
elseif sol(3) < 1e-2 && sol(2) > 1e-2
    sol(3) =  0.01;
end


fittedParams    = sol;
fittedQuantiles = two_part_normal_inverse_cdf(quantiles, sol(1), sol(2), sol(3));

end