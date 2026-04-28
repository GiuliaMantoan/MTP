function tpnVariance = two_part_normal_variance(parameters)

sigma1 =  parameters(:,2);
sigma2 =  parameters(:,3);

tpnVariance = (1-2/pi)*(sigma2-sigma1).^2+sigma1.*sigma2;

end