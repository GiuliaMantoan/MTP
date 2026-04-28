function result = gauss_legendre_quadrature(f, a, b, nodes)
    % Define the nodes and weights for Gauss-Legendre quadrature
    switch nodes
        case 2
            n = [-1/sqrt(3) 1/sqrt(3)]; % location values for n=2
            w = [1 1]; % weight for n=2
        case 4
            n = [-0.86113 -0.33998 0.33998 0.86113]; % location values for n=4
            w = [0.34785 0.625214 0.625214 0.34785]; % weights for n=4
        case 6
            n = [-0.93246 -0.66120 -0.23861 0.23861 0.66120 0.93246]; % location values for n=6
            w = [0.17132 0.36076 0.46791 0.46791 0.36076 0.17132]; % weights for n=6
        otherwise
            error('Supported node values are 2, 4, and 6.');
    end

    % Map the nodes from the interval [-1, 1] to [a, b]
    x = ((b - a) / 2) * n + ((b + a) / 2);
    % Adjust the weights for the interval [a, b]
    adjusted_weights = (b - a) / 2 * w;
    % Compute the integral using Gauss-Legendre quadrature
    result = sum(adjusted_weights .* f(x));
end


