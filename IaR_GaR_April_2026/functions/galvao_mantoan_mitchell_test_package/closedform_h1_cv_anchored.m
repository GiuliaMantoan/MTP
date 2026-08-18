function [Kv_sim_h1, CVMv_sim_h1] = closedform_h1_cv_anchored(rvec, N, seed, KS_tab, CVM_tab)
% 2026-08-13: same as closedform_h1_cv.m (simulates the limiting
% Brownian-bridge distribution for horizon 1) but RESCALES the resulting
% draws so their OWN 95th percentile lands exactly on R&S(2019)'s
% published tabulated constants (KS_tab, CVM_tab), instead of whatever
% our own simulation happens to produce. Our self-simulated critical
% value depends on rvec grid resolution and N (multipleH_smallH_P500_
% diagnostic_h1rankfix.m found ours coming in at ~1.30-1.32 vs R&S's
% published 1.339, likely grid-resolution noise) -- rescaling anchors the
% marginal test to hit the literal tabulated number exactly, while still
% producing a full distribution of draws for
% CVfinalbootstrsuptest_shifted2_h1rankfix.m's rank-preserving quantile
% mapping to map onto (so the Sup Test's own joint critical value is
% corrected too, not just the marginal readout -- a single constant alone
% can't do that, see that function's header).

if nargin < 4 || isempty(KS_tab)
    KS_tab = 1.339;
end
if nargin < 5 || isempty(CVM_tab)
    CVM_tab = 0.460;
end

[Kv_sim_h1, CVMv_sim_h1] = closedform_h1_cv(rvec, N, seed);

N_sim = numel(Kv_sim_h1);
scaleK = KS_tab / Kv_sim_h1(ceil(N_sim*0.95));
scaleM = CVM_tab / CVMv_sim_h1(ceil(N_sim*0.95));

Kv_sim_h1 = Kv_sim_h1 * scaleK;
CVMv_sim_h1 = CVMv_sim_h1 * scaleM;
end
