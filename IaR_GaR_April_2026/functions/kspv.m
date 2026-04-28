function z = kspv(d, n)

    % Calculate lambda
    lambda = (sqrt(n) + 0.12 + (0.11 / sqrt(n))) * d;
    
    % Generate sequence j
    j = 1:100;
    
    % Calculate z
    z = 2 * sum(((-1).^(j-1)) .* exp(-2 * (j.^2) * (lambda^2)));
    
end