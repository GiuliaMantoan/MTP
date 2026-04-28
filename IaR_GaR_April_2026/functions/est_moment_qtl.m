function [mu, sigma, skw, kur] = est_moment_qtl(quantiles)

    % lenght check
    if size(quantiles,2) == 1
        error('Transpose the vector of quantile');
    end
    
    % Define the weight 
    weight = 1 / size(quantiles,2);
    
    % Compute mu: weighted sum of quantiles
    mu = sum(quantiles .* weight);
    
    % Compute sigma^2: weighted sum of squared differences from the mean
    sigma_2 = sum(((quantiles - mu).^2) .* weight);
    sigma = sqrt(sigma_2);  % sqrt of sigma^2
    
    % Compute skw: weighted sum of (quantiles - mu)^3 / sigma
    skw = sum((((quantiles - mu) / sigma).^3) .* weight);

    % Compute kurtosis: weighted sum of (quantiles - mu)^4 / sigma^4
    kur = sum((((quantiles - mu) / sigma).^4) .* weight);

end