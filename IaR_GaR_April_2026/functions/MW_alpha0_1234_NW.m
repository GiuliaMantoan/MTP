function [stat_uncens] = MW_alpha0_1234_NW(x, lags, prewhite, z_l, z_u)

    % The code is adapted from Knuppel 2015 symmetric test. The PITs are taken
    % as given in this case and denoted by x. To symmetrise the PITs,
    % standardised-PIT random variables are constructed with same set of
    % moments as before. 
    % A separate test statistic for the frequency of outtakes outside the
    % censored region is also added. 


    % Moment test for time series data, using Newey-West covariance matrix
    % restricted to the censored region [z_l, z_u].
    if nargin < 4
       z_l = 0.05;
       z_u = 0.95;
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

    % Identify which observations lie in the interval [z_l, z_u].
    idx_in = (x >= z_l & x <= z_u);

    % Only use x values in [z_l, z_u] to compute s_pit and the moment conditions
    x_in  = x(idx_in);
    N_in  = sum(idx_in);  % number of in-range observations

    if N_in == 0
       error('No x values lie in the interval [z_l, z_u]. Cannot compute test.');
    end
    
    % Compute symmetrized PIT only on the in-range x_in
    s_pit_in = sqrt(12 / (z_u - z_l)^2) .* (x_in - 0.5 * (z_l + z_u));
    
    % z has four columns (1st, 2nd, 3rd, 4th "moments")
    z_in = [ s_pit_in, ...
             (s_pit_in.^2 - 1), ...
             (s_pit_in.^3), ...
             (s_pit_in.^4 - 1.8) ];

    %-----------------------------
    % Compute odd-moment part
    %-----------------------------
    z_odd        = z_in(:, [1 3]);         % columns for the 1st and 3rd moments
    z_odd_mean   = sum(z_odd, 1) .* (N_in^(-0.5));  % scaled sum
    [phi_odd, ~] = nwX(z_odd, prewhite, lags);     % NW covariance on the full Nx2 matrix
    stat_odd     = z_odd_mean * (phi_odd \ z_odd_mean'); 
    
    %-----------------------------
    % Compute even-moment part
    %-----------------------------
    z_even        = z_in(:, [2 4]);        % columns for the 2nd and 4th moments
    z_even_mean   = sum(z_even, 1) .* (N_in^(-0.5)); 
    [phi_even, ~] = nwX(z_even, prewhite, lags); 
    stat_even     = z_even_mean * (phi_even \ z_even_mean');  
    
    %-----------------------------
    % Final statistic
    %-----------------------------
    stat_uncens = stat_odd + stat_even;
end
