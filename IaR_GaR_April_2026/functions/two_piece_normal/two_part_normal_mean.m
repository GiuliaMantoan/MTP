function tpnMean = two_part_normal_mean(parameters)

mode = parameters(:,1);
% [sigma1,sigma2] = reparametrise_two_part_normal(parameters(:,2),parameters(:,3));
sigma1 =  parameters(:,2);
sigma2 =  parameters(:,3);

tpnMean = mode+sqrt(2/pi)*(sigma2-sigma1);
% tpnMean(parameters(:,2) == 0) = parameters(parameters(:,2) == 0,1);
tpnMean(parameters(:,2)<1e-4 | parameters(:,3)<1e-4) = parameters(parameters(:,2)<1e-4 | parameters(:,3)<1e-4, 1);

end