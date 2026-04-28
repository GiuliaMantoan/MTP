function [fittedParams,fittedQuantiles] = fit_tpn_to_quantiles(data,quantiles,weights)

% INPUT:
% data: is a T x M vector where T is the number of time periods and M the
%      number of quantiles
% quantiles: is 1 x M vector of quantiles. For example: [0.05 ... 0.5 ...
%           0.95]
% weights: weights provided to a given quantile in the minimization pb 

% OUTPT
% fittedParams: parameters from the 2-piece normal 
% fittedQuantiles: quantiles fitted from the minimization pb 

fittedParams = zeros(size(data,1),3);
fittedQuantiles = zeros(size(data));
paramsInit = [0 1 1]; % initial guesses

for tt = 1:size(data,1)
    
    distToFit = @(x) (two_part_normal_inverse_cdf(quantiles,x(1),x(2),x(3))-data(tt,:)).*weights;
    temp = fsolve(distToFit,paramsInit,optimset('MaxFunEval',20000,'MaxIter',20000,'TolFun',1e-6,'TolX',1e-20));
    
    % modified by SM June 2025.
    % in the original version of the code avaiable in Simon's website
    % colomn 2 and 3 of fittedParams where filled with the squre root of
    % dispersion and the skewness in BOE formulation (what in wiki is label
    % as xi

    fittedParams(tt,1) = temp(1);  %mode
    %fittedParams(tt,2) = sqrt(2*(temp(2)^2)*(temp(3)^2)/(temp(2)^2+temp(3)^2));  %unc
    %fittedParams(tt,3) = sqrt(2/pi)*(temp(3)-temp(2)); %skew
    fittedParams(tt,2) = temp(2);
    fittedParams(tt,3) = temp(3);    

    fittedQuantiles(tt,:) = two_part_normal_inverse_cdf(quantiles,...
        fittedParams(tt,1),temp(2),temp(3));
%     paramsInit = temp;

    
end

end