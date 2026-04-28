function logL = ar1logli(y, mu, rho, sigma)

    % Number of observations
    n = length(y);

    % First term
    term1 = -0.5 * log(2 * pi);

    % Second term
    term2 = -0.5 * log(sigma / (1 - rho^2));

    % Third term
    term3 = -((y(1) - (mu / (1 - rho)))^2) / (2 * sigma / (1 - rho^2));

    % Fourth term
    term4 = -((n - 1) / 2) * log(2 * pi);

    % Fifth term
    term5 = -((n - 1) / 2) * log(sigma);

    % Sixth term
    residuals = y(2:n) - mu - rho * y(1:n-1);
    term6 = -sum((residuals.^2) / (2 * sigma));

    % Total log-likelihood
    logL = term1 + term2 + term3 + term4 + term5 + term6;
    
end
