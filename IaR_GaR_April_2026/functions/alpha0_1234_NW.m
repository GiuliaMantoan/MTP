function [stat,pval] = alpha0_1234_NW(x,lags,prewhite)
    
    % remove NaN
    x = x(~isnan(x));
    x = x';
    
    % Moment test for time series data, using Newey-West covariance matrix

    if nargin < 3
        prewhite = 0;
        if nargin < 2
            lags = -1;
        end
    end
    % x are the PITs here
    s_pit = sqrt(12) * (x - 0.5);
    z = [s_pit (s_pit.^2 - 1) (s_pit.^3) (s_pit.^4 - 1.8)];
    z_orig = z;

    z = z_orig(:, [1 3]);
    y = sum(z) * (size(z,1)^(-0.5));
    [phi,~] = nwX(z, prewhite, lags);       
    stat_odd = y * (phi\ y');

    z = z_orig(:, [2 4]);
    y = sum(z) * (size(z,1)^(-0.5));
    [phi,~] = nwX(z, prewhite, lags);       
    stat_even = y *(phi \ y');

    df = 4;
    stat = stat_odd + stat_even;
    pval = 1 - chi2cdf(stat, df);
end
