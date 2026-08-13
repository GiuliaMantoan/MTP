function [Kv_sim, CVMv_sim] = closedform_h1_cv(rvec, N, seed)
% 2026-08-13: simulates the LIMITING (Brownian bridge) distribution of
% the unweighted KS/CvM statistics for horizon 1, matching R&S(2019)'s
% own CVfinalALLhac.m (found in their published replication package,
% RossiSekhposyanJoE2018Codes/MonteCarlos/) -- NOT a block bootstrap.
%
% Horizon 1's PIT is genuinely iid under the null (u_t^(1)=eps_{t+1}, no
% MA terms), so its limiting empirical process is a plain Brownian bridge
% with covariance kernel Sigma(r1,r2) = min(r1,r2) - r1*r2 -- no nuisance
% parameters, no need to estimate anything from data. R&S(2019) use this
% closed-form/simulated-once distribution for h=1 instead of bootstrapping
% (see their Table 1 and Section 3). This function reproduces that
% simulation so callers can extract whatever quantile they need (e.g.
% the 95th percentile for a standalone test, or the 1-alpha/H percentile
% for a Bonferroni-corrected test at a given H) rather than hard-coding a
% single critical value.
%
% Returns N draws of Kv (=max|v(r)|) and CVMv (=mean(v(r)^2)) from this
% limiting process, sorted ascending, ready for quantile extraction via
% e.g. Kv_sim(ceil(N*0.95)).

if nargin < 3
    stream1 = RandStream('mrg32k3a','seed',987654);
else
    stream1 = RandStream('mrg32k3a','seed',seed);
end

[rx, ry] = meshgrid(rvec, rvec);
Sigma = min(rx,ry) - rx.*ry;
L = cholcov(Sigma);   % L' * L = Sigma (cholcov returns L such that L'*L=Sigma, matching R&S(2019)'s own usage)

draws = (L' * randn(stream1, size(L,1), N));   % size(rvec,2) x N
Qvabs = abs(draws);
Kv_sim = sort(max(Qvabs,[],1), 'ascend')';
CVMv_sim = sort(mean(draws.^2,1), 'ascend')';
end
