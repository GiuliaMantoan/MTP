function [fittedParams,fittedQuantiles] = fit_tpn_to_quantiles_SM_reparam(data,quantiles,mode0, sigma1_0, sigma2_0)

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
[~,jq05] = min(abs(quantiles-.05)); % lower qt
[~,jq95] = min(abs(quantiles-.95)); % upper qt

% Initial guess based on empirical quantiles
if nargin < 3

    mode0    = data(jq50); % mode
    
    z95      = norminv(0.95, 0, 1);
    z05      = norminv(0.05, 0, 1);

    sigma1_0 = (data(jq50) - data(jq05)) / abs(z05); % s1
    sigma2_0 = (data(jq95) - data(jq50)) / z95; %s2

end

mode0 = 1;
lambda1_0 = 0;
lambda2_0 = 0;


paramsInit = [mode0, lambda1_0, lambda2_0]; % initial guess

% Bounds: ensure sigmas positive
LB = [-Inf, -Inf , -Inf];
UB = [ Inf, Inf, Inf];

opts = optimoptions('fminunc', ...
                        'Algorithm','quasi-newton', ...
                        'Display','off', ...
                        'MaxFunctionEvaluations',1e4, ...
                        'MaxIterations',1e4, ...
                        'OptimalityTolerance',1e-6, ...
                        'StepTolerance',1e-6);


resid = @(x) ( two_part_normal_inverse_cdf_reparam(quantiles([jq05, jq50, jq95]), x(1), x(2), x(3)) - data([jq05, jq50, jq95]) ); % residual from actual qt and qt obtained with the parameters
objfun = @(x) norm(resid(x)); % L2 norm 

sol = fminunc(objfun, paramsInit, opts);

ss1 = log(sol(2));
ss2 = log(sol(3));

fittedParams    = [sol(1), ss1, ss2];
fittedQuantiles = two_part_normal_inverse_cdf(quantiles, sol(1), ss1, ss2);



end