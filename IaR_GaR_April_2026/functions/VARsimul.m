function [XXb, Sigma] = VARsimul(X, p, nSim)
%% VARsimul: fit VAR to data and generate samples of equal length
%
% Description: VARsimul fits a VAR with p lags to the provided data X and
% uses the estimated parameters to simulate samples of equal length (using
% the same p initial values contained in X and Gaussian shocks with
% covariance matrix estimated from the data).
% 
% Input arguments:
% - X : T-by-N matrix containing data used to estimate the VAR parameters.
% - p : Number of lags to include in the estimated VAR.
% - nSim : Number of samples to generate using the estimated parameters.
%
% Output arguments:
% - XXb : T-by-N-by-nSim array, such that XXb(:,:,m) contains a simulated
%         sample.




%% Estimate VAR parameters using data provided in X
[T, N] = size(X);
% Construct Z matrix of right-hand side variables (constant plus lags of X)
Z = ones(T - p, 1);
for j = 1:p
    Z = [Z, X((p + 1 - j):(T - j), :)];
end
y = X((p + 1):end, :);
% Estimate VAR coefficient matrix and innovation covariance matrix
beta = Z \ y;
Sigma = cov(y - Z * beta);

%% Simulate samples
XXb = zeros(T, N, nSim); % array to store simulated samples
for m = 1:nSim
    clear('Xb')
    Xb(1:p, :) = X(1:p, :); % use the p initial values from the data X
    err = randn(T, N) * chol(Sigma); % draw errors
    for t = (p + 1):T
        % Construct Zt matrix of time t values for conditioning variables
        Zt = 1;
        for j = 1:p
            Zt = [Zt, Xb(t - j, :)];
        end
        Xb(t, :) = Zt * beta + err(t, :);
    end
    XXb(:, :, m) = Xb;
end