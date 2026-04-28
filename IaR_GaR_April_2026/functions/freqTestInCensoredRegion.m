function [stat_freq] = freqTestInCensoredRegion(x, lags, prewhite, z_l, z_u)
% Test whether the fraction of x in [z_l, z_u] equals its nominal coverage alpha=(z_u-z_l).
% * D_t = 1{x_t in [z_l,z_u]} - alpha
% * We average D_t, do a HAC-based standard error, and form a chi-square(1) test.

    if nargin < 5
        z_u = 0.95;
        z_l = 0.05;
    end
    if nargin < 3
        prewhite = 0;
        if nargin < 2
            lags = -1;
        end
    end
 
    % remove NaN
    x = x(~isnan(x));
    x = x';
    
    % Nominal best coverage region size
    alpha = z_u - z_l;
    
    % Construct D
    N = length(x);
    idx_in = (x >= z_l & x <= z_u);
    D = double(idx_in) - alpha;

    D_bar = sum(D)/sqrt(N);

    % HAC estimate of Var(D) which is univariate here
    [phiD, ~] = nwX(D, prewhite, lags);

    stat_freq = (D_bar^2) / phiD;
end
