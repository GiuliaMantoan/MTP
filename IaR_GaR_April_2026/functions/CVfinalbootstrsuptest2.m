function [result_sup, result_bonf, result_h] = CVfinalbootstrsuptest2(el, bootMC, pit, rvec, w, wi, stream1)
% CVfinalbootstrsuptest2
% Computes bootstrap critical values for:
%   - Sup-type KS and CVM tests (unweighted, weighted by w, weighted by wi)\%   - Bonferroni-adjusted horizon-by-horizon tests
%   - Horizon-specific critical values
%
% Inputs:
%   el     : block length for stationary bootstrap
%   bootMC : number of bootstrap draws
%   pit    : [P x H] matrix of PIT values (z-statistics) for each time and horizon
%   rvec   : vector of thresholds in [0,1] for empirical CDF comparisons
%   w      : 1xH vector of sup-test weights
%   wi     : 1xH vector of alternative sup-test weights
%   stream1: RandStream object controlling random numbers
%
% Outputs:
%   result_sup  : 3x2 matrix of 95%% critical values for sup-type tests
%                 rows correspond to {max, w, wi}; columns to {KS, CVM}
%   result_bonf : 2xH matrix of Bonferroni-adjusted critical values
%                 row 1: KS, row 2: CVM, columns for each horizon
%   result_h    : 2xH matrix of horizon-specific critical values
%                 row 1: KS, row 2: CVM, columns for each horizon

% Number of horizons
H = size(pit,2);

% Preallocate storage for bootstrap statistics
KSv          = zeros(bootMC,H);    % KS stat per bootstrap and horizon
CVMv         = zeros(bootMC,H);    % CVM stat per bootstrap and horizon
KSv_sup_max  = zeros(bootMC,1);    % unweighted sup KS
CVMv_sup_max = zeros(bootMC,1);    % unweighted sup CVM
KSv_sup_w    = zeros(bootMC,1);    % weighted w KS
CVMv_sup_w   = zeros(bootMC,1);    % weighted w CVM
KSv_sup_wi   = zeros(bootMC,1);    % weighted wi KS
CVMv_sup_wi  = zeros(bootMC,1);    % weighted wi CVM

P        = size(pit,1);            % total out-of-sample observations
nThresh  = numel(rvec);            % number of thresholds

% Main bootstrap loop
tic;

for b = 1:bootMC
    % Temporary bootstrap sample of K* statistics per horizon x threshold
    K_star = zeros(H, nThresh);

    for h = 1:H
        % Generate i.i.d. multipliers for stationary bootstrap increment
        z = (1/sqrt(el)) * randn(stream1, P-el+1, 1);

        % Empirical CDF differences:
        % indicator whether PIT <= thresholds, minus mean CDF
        emp_cdf      = (repmat(pit(:,h),1,nThresh) <= repmat(rvec,P,1));
        mean_emp_cdf = mean(emp_cdf,1);

        % Accumulate block sums scaled by multiplier z
        for j = 1:(P-el+1)
            block_diff        = emp_cdf(j:j+el-1,:) - repmat(mean_emp_cdf,el,1);
            K_star(h,:)       = K_star(h,:) + (P^(-1/2) * z(j)) .* sum(block_diff,1);
        end

        % KS: max absolute deviation; CVM: average squared deviation
        KSv(b,h)   = max(abs(K_star(h,:)));
        CVMv(b,h)  = mean(K_star(h,:).^2);
    end

    % Sup-type statistics across horizons
    KSv_sup_max(b)  = max(KSv(b,:));
    CVMv_sup_max(b) = max(CVMv(b,:));
    KSv_sup_w(b)    = max(KSv(b,:) .* w);
    CVMv_sup_w(b)   = max(CVMv(b,:) .* w);
    KSv_sup_wi(b)   = max(KSv(b,:) .* wi);
    CVMv_sup_wi(b)  = max(CVMv(b,:) .* wi);
end

timeElapsed = toc; %#ok<NASGU>

% Sort bootstrap distributions
KSv          = sort(KSv,1);
CVMv         = sort(CVMv,1);
KSv_sup_max  = sort(KSv_sup_max);
CVMv_sup_max = sort(CVMv_sup_max);
KSv_sup_w    = sort(KSv_sup_w);
CVMv_sup_w   = sort(CVMv_sup_w);
KSv_sup_wi   = sort(KSv_sup_wi);
CVMv_sup_wi  = sort(CVMv_sup_wi);

% Compute 95%% critical values
alpha = 0.05;
cvKv_h    = KSv(round(bootMC*(1-alpha),1),:);      % horizon-specific KS
cvMv_h    = CVMv(round(bootMC*(1-alpha),1),:);     % horizon-specific CVM
cvKv_max  = KSv_sup_max(round(bootMC*(1-alpha),1));
cvMv_max  = CVMv_sup_max(round(bootMC*(1-alpha),1));
cvKv_w    = KSv_sup_w(round(bootMC*(1-alpha),1));
cvMv_w    = CVMv_sup_w(round(bootMC*(1-alpha),1));
cvKv_wi   = KSv_sup_wi(round(bootMC*(1-alpha),1));
cvMv_wi   = CVMv_sup_wi(round(bootMC*(1-alpha),1));

% Bonferroni-adjusted thresholds per horizon
bonf_idx = round(bootMC * (1 - alpha/H));
cvKv_bonf = KSv(bonf_idx,:);
cvMv_bonf = CVMv(bonf_idx,:);

% Package outputs
result_sup  = [cvKv_max, cvMv_max; cvKv_w, cvMv_w; cvKv_wi, cvMv_wi];
result_bonf = [cvKv_bonf; cvMv_bonf];
result_h    = [cvKv_h; cvMv_h];
end
