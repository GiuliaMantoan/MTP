function [x, w] = lgwt(N, a, b)
    % lgwt computes the nodes x and weights w for Gauss-Legendre quadrature
    % over the interval [a, b] with N points.
    
    % Calculation of the zeros of the Legendre polynomial
    x = cos(pi * (4 * (1:N) - 1) ./ (4 * N + 2))';
    P0 = zeros(N, 1);
    P1 = ones(N, 1);
    
    for k = 2:N
        P2 = ((2 * k - 1) * x .* P1 - (k - 1) * P0) / k;
        P0 = P1;
        P1 = P2;
    end
    
    x = x - P1 ./ (2 * P2);
    
    % Linear map from [-1,1] to [a,b]
    x = (a * (1 - x) + b * (1 + x)) / 2;
    
    % Compute the weights
    w = (b - a) ./ ((1 - x.^2) .* (polyval(polyder(legendrePoly(N)), x)).^2) * (2 / N)^2;
end