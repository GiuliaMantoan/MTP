function [test1, test2, pval1, pval2] = berk(z)

    % Remove missing values from z 
    z = z(~isnan(z));
    z = z';


    % Apply inverse CDF of the normal distribution to each element in z 
    for i = 1:size(z, 1)
        for j = 1:size(z, 2)
            z(i, j) = norminv(z(i, j));
        end
    end
    
lagged_z = lagmatrix(z, 1); % Create a lagged version of z with lag of 1

x = [ones(size(z, 1) , 1), lagged_z]; % Combine the column of ones with the lagged z

% Drop the first row of z and x
y = z(2:end, :); % Drop the first row of z
x = x(2:end, :); % Drop the first row of x

b = x \ y; % Equivalent to b = (x' * x) \ (x' * y)

% Calculate residuals
residuals = y - x * b;

% Compute the variance-covariance matrix of the residuals
sigma = cov(residuals);


mu = b(1); % coeff constant 
rho = b(2); % coeff AR

logL_full = ar1logli(z, mu, rho, sigma); % Log-likelihood with mu=mu, rho=rho, sigma=sigma
logL0 = ar1logli(z, 0, 0, 1);  % Log-likelihood with mu=0, rho=0, sigma=1
logL1 = ar1logli(z, mu, 0, sigma); % Log-likelihood with mu=mu, rho=0, sigma=sigma

% Calculate the test statistics
test1 = -2 * (logL0 - logL_full);
test2 = -2 * (logL1 - logL_full);

% If test1 is complex, replace both with NaN
if ~isreal(test1)
    test1 = NaN;
    test2 = NaN;
end

% Calculate the p-values using chi2cdf
pval1 = chi2cdf(test1, 3, 'upper'); 
pval2 = chi2cdf(test2, 1, 'upper');


end






