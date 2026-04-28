function inverseDensity = two_part_normal_inverse_cdf(F,mode,sigma1,...
    sigma2)
% This function returns the inverse cumulative density function for a two 
% part normal density
% This function takes as inputs the parameters for a given two part normal
% distribution (in the form of mode and the two different standard
% deviations, one below and one above the mode), and the probability F such that
% the value which has F probability mass below it is returned. This function outputs a
% single value.

% Self-explanatory error as probability value b
if any(F<0) || any(F>1)
    error(['Inputs to two_part_normal_inverse_cdf must ',...
        'be valid cumulative densities (between 0 and 1).']);
end   

% The calculations below exploit the fact that the two-part normal
% distribution is basically a normal distribution below the mode and a
% different normal distribution below the mode. 
% So say that 55% of the
% probability mass lies below the mode (this is the value of the variable
% 'cutoff' defined below). Then there are different procedures to compute
% the value depending on whether F is greater or less than 0.55, but they
% both exploit the fact that a normal cdf can be used together with the
% parameters for this distribution.

cutoff         = sigma1/(sigma1+sigma2);
lessThanCutoff = (F<=cutoff);
lowF           = F.*lessThanCutoff + cutoff*(~lessThanCutoff);
highF          = F.*(~lessThanCutoff) + cutoff*lessThanCutoff;
 
inverseDensity = ...
    lessThanCutoff.*norminv((sigma1+sigma2)/(2*sigma1)*lowF,mode,sigma1) + ...
    ~lessThanCutoff.*norminv((highF-(sigma1-sigma2)/(sigma1+sigma2))*...
    (sigma1+sigma2)/(2*sigma2),mode,sigma2);