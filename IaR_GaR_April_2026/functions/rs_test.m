function [rs_boot_crit_value, rs_test_stats, rs_logic, KS, CVM] = rs_test(z, rvec)
% INPUT
% z - T by 1 vector of the probability integral transforms
% rvec - 1 by x vector with values between 0 and 1 guiding the
% discretization of the uniform
% *************************************************************************
% Implementing the test
% *************************************************************************
% Creating the test statistics
% cumcumz corresponds to the xi and phi on page 8, respectively
T = size(z,1);
cumcumz = (repmat(z,1,size(rvec,2)) < repmat(rvec,size(z,1),1)) - repmat(rvec,size(z,1),1);
v = sum(cumcumz,1)/sqrt(T);

KS = max(abs(v)); % Kolmogorov-Smirnov test statistics
CVM = mean(v.^2); % Cramer-von-Mises test statistics

rs_test_stats.ks_stat  = KS;
rs_test_stats.cvm_stat = CVM;

% Simulate the critical values applicable for multi-step-ahead density
% forecast calibration
% Implements the test based on the Theorem 3
% setup some preliminaries for the bootstrap
el = floor(T.^(1/3)); % block length
bootMC = 200;

tableboot = CVfinalbootstrapInoue(el, bootMC, z, rvec);
rs_boot_crit_value.ks = tableboot(:,1)';
rs_boot_crit_value.cvm = tableboot(:,2)';
KS_logic = (KS > rs_boot_crit_value.ks); 
CVM_logic = (CVM > rs_boot_crit_value.cvm);

rs_logic.ks_array = KS_logic;
rs_logic.cvm_array = CVM_logic;

% *************************************************************************
% Report the results
% *************************************************************************
delete Results.out; diary Results.out;
% disp('Test statistics [KS CVM]');
% disp([KS CVM]);
% disp('Critical Values, 1% 5% 10%')
% disp('Multi-step-ahead density calibration test');
% disp('Kolmogorov-Smirnov Test');
% disp(rs_boot_crit_value.ks);
% disp('Cramer-von-Mises Test');
% disp(rs_boot_crit_value.cvm);
% % Decision logic for KS test
% if KS < rs_boot_crit_value.ks(3)
%     disp('KS Test: Fail to reject the null at 10%, 5%, and 1%.');
% elseif KS < rs_boot_crit_value.ks(2)
%     disp('KS Test: Reject the null at 10%, but fail to reject at 5% and 1%.');
% elseif KS < rs_boot_crit_value.ks(1)
%     disp('KS Test: Reject the null at 10% and 5%, but fail to reject at 1%.');
% else
%     disp('KS Test: Reject the null at all significance levels (1%, 5%, and 10%).');
% end
% 
% % Decision logic for CVM test
% if CVM < rs_boot_crit_value.cvm(3)
%     disp('CVM Test: Fail to reject the null at 10%, 5%, and 1%.');
% elseif CVM < rs_boot_crit_value.cvm(2)
%     disp('CVM Test: Reject the null at 10%, but fail to reject at 5% and 1%.');
% elseif CVM < rs_boot_crit_value.cvm(1)
%     disp('CVM Test: Reject the null at 10% and 5%, but fail to reject at 1%.');
% else
%     disp('CVM Test: Reject the null at all significance levels (1%, 5%, and 10%).');
% end

diary off;
end
