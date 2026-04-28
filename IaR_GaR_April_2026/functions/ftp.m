function yy = ftp(y, m, s1, s2)
    
% Initialize yy as a copy of y
yy = y;
    
% Compute the scaling factor
fac = ((pi * ((s1 + s2)^2)) / 2)^(-0.5);
    
 % Iterate over rows and columns of y
 for i = 1:size(y, 1)
     for j = 1:size(y, 2)
            if y(i, j) < m
                yy(i, j) = fac * exp((-0.5 * ((y(i, j) - m)^2)) / (s1^2));
            elseif y(i, j) == m
                yy(i, j) = fac;
            else
                yy(i, j) = fac * exp((-0.5 * ((y(i, j) - m)^2)) / (s2^2));
            end
      end
  end
end
