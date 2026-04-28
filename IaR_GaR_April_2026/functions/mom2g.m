function sig = mom2g(s, gam)

    % Compute alpha and beta
    alpha = sqrt(1 - gam);
    beta = sqrt(1 + gam);
    
    % Compute sig
    sig = ((1 - (2 / pi)) * (((s / alpha) - (s / beta))^2)) + ((s^2) / (beta * alpha));
    
end