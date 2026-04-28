function y = rsn(n, loc, scale, shape)
    % Generate n random numbers from a skew-normal distribution
    % shape: Skewness parameter (alpha)
    % loc: Location parameter (mu)
    % scale: Scale parameter (sigma)
    
    % Generate standard normal random variables
    u0 = randn(n, 1);
    v = randn(n, 1);
    
    % Apply the skew transformation
    delta = shape / sqrt(1 + shape^2);
    u1 = delta * abs(u0) + sqrt(1 - delta^2) * v;
    
    % Transform to skew-normal distribution
    y = loc + scale * u1;
    y = y';
end