function P = legendrePoly(n)
    % legendrePoly computes the Legendre polynomial of degree n
    P = zeros(1, n + 1);
    P(end) = 1;
    
    for k = n:-1:1
        P = [0 P] - [(2*k-1)/(k) P(1:end-1)] * 2;
    end
    
    P = P / P(end);
end