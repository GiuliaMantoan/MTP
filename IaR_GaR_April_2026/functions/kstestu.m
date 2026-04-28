function [d, p_value] = kstestu(x)


d = max(abs(mean(bsxfun(@le, x, x'))' - x)); % comparison of the empirical cdf with the theoretical cdf

p_value =  kspv(d, size(x,1));

end