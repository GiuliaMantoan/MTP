function result = gaussLegendreQuadrature(func, a, b, n)
    % func - the function to be integrated
    % a, b - the limits of integration
    % n - the number of points for Gauss-Legendre quadrature
    
    % Obtain the nodes and weights for the Gauss-Legendre quadrature
    [x, w] = lgwt(n, a, b);
    
    % Compute the integral
    result = sum(w .* func(x));
end