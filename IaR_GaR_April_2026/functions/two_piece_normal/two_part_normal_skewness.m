function tpnSkew = two_part_normal_skewness(parameters)

sigma1 =  parameters(:,2);
sigma2 =  parameters(:,3);

tpnSkew = sqrt(2/pi) * (sigma2 - sigma1) * ( (4/pi - 1) * (sigma2 - sigma1)^2 + sigma1 * sigma2 );

end