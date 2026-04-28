function [parm, sig1, sig2] = momtopar(m, mu, v)
    % Initial assignments
    parm = mu;
    alpha = (m - mu) * sqrt(pi / 2);
    
    % Calculate discriminant
    disc = sqrt(((1 - 4 * (1 - (2 / pi))) * (alpha^2)) + 4 * v);
    
    % Conditional assignment for sig1
    if (disc - alpha) > 0
        sig1 = (disc - alpha) / 2;
    else
        sig1 = -1 * (disc + alpha) / 2;
    end
    
    % Calculate sig2
    sig2 = alpha + sig1;
end