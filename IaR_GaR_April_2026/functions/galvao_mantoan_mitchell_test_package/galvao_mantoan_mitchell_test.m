%% Galvao-Mantoan-Mitchell -- Sup Test / Bonferroni, applied to real (fan chart) data
%
% 2026-08-13: rewritten to use the latest bootstrap developed in the
% Path Evaluation Monte Carlo project, in place of the original
% CVfinalbootstrsuptest2 (which modeled every horizon's bootstrap column
% independently, with no cross-horizon correlation at all).
%
% Two fixes are ported over from that project's validated Monte Carlo work:
%
% 1. CROSS-HORIZON CORRELATION: horizons that share the same underlying
%    shock are index-shift aligned in the bootstrap (horizon h's block i
%    reads a shared eta sequence at position i+h, so two horizons that
%    share a dominant shock also share the same eta draw for it). This
%    project's Monte Carlo work also derived an ANALYTICALLY-WEIGHTED
%    version of this (correlation matching the true DGP's theta-driven
%    decay exactly) -- but that weighting is derived FROM the simulated
%    IMA(1) DGP's known theta parameter, which has no meaningful value
%    for real data with an unknown true dependence structure. This script
%    therefore uses the FULL-SHARING version (correlation=1 at the
%    aligned shift, no theta needed) -- a defensible, model-free choice:
%    it doesn't assume zero correlation like the original, and doesn't
%    assume a specific wrong decay parameter either.
%
% 2. HORIZON 1: h=1's PIT has no MA dependence under the null of correct
%    specification (true for ANY correctly-specified model, not just the
%    simulated DGP) -- R&S(2019) never bootstrap it, using the closed-form
%    Brownian-bridge asymptotic distribution instead. Applying the block
%    bootstrap to horizon 1 anyway was found to bias its critical value
%    (both the marginal test AND, through the joint max, the Sup Test's
%    own critical value). Fix: horizon 1's raw bootstrap column is
%    rank-preserving quantile-mapped onto the closed-form distribution,
%    anchored so its own 95th percentile matches R&S(2019)'s published
%    tabulated constants (KS=1.339, CvM=0.460) exactly.
%
% FILES NEEDED alongside this script (self-contained, no other project
% dependencies):
%   closedform_h1_cv.m
%   closedform_h1_cv_anchored.m
%   CVfinalbootstrsuptest_shifted_h1rankfix.m
%   size_statistic_h2_shifted_h1rankfix.m
%
% Everything else (nModels, models, zinf_all, results.(mod).KS_vec /
% CVM_vec already populated from the observed data) is assumed to already
% exist in the workspace exactly as in the original script -- this only
% replaces the bootstrap step and its output fields
% (results.(mod).gmm_ks / gmm_cvm / gmm_ks_bonf / gmm_cvm_bonf), so
% nothing downstream needs to change.

% 2026-08-13: MC/bootMC updated to match the Path Evaluation project's
% current production defaults (multipleH_size.m) -- bootMC=200 was found
% to give the same size as bootMC=5000 (multipleH_size_bootMC_diagnostic.m),
% so 200 is used for speed; MC=5000 for a tight read on the resulting
% rejection frequency.
MC     = 5000;
bootMC = 200;
rng(bootMC, 'twister');
rvec_gmm = 0:0.001:1;

% ---- Precompute the closed-form horizon-1 distribution ONCE -- it only
% depends on rvec_gmm, not on the model or the data, so this must NOT be
% recomputed inside the loop over models or the bootstrap loop. Anchored
% to R&S(2019)'s own published 5% critical values (KS=1.339, CvM=0.460).
N_sim_h1 = 200000;
fprintf('Precomputing horizon-1 closed-form (Brownian bridge) distribution,\n');
fprintf('anchored to R&S(2019)''s tabulated KS=1.339 / CvM=0.460...\n');
[Kv_sim_h1, CVMv_sim_h1] = closedform_h1_cv_anchored(rvec_gmm, N_sim_h1, 987654);
fprintf('Done.\n\n');

for m = 1:nModels
    mod  = models{m};
    zinf = zinf_all{m};
    z    = zinf';
    P    = size(z, 1);
    % 2026-08-13: switched from P^(1/4) to P^(1/3) -- R&S(2019)'s own
    % baseline block length (P^(1/4) is only their robustness-panel
    % choice). The Path Evaluation project's own Monte Carlo work found
    % P^(1/3) gives Sup Test size consistently closer to nominal 5% than
    % P^(1/4), matching the choice now used in that project's production
    % sweep (multipleH_size.m).
    el   = floor(P^(1/3));
    Hz   = size(z, 2);
    KS   = results.(mod).KS_vec(:)';
    CVM  = results.(mod).CVM_vec(:)';
    QVrejvecs       = zeros(MC, 3);
    CVMrejvecs      = zeros(MC, 3);
    QVrejvecs_bonf  = zeros(MC, 3);
    CVMrejvecs_bonf = zeros(MC, 3);
    parfor j = 1:MC
        stream1 = RandStream('mrg32k3a', 'seed', 4829575);
        stream1.Substream = j;
        [QVrej_j, CVMrej_j, QVbonf_j, CVMbonf_j] = ...
            size_statistic_h2_shifted_h1rankfix(z, KS, CVM, Hz, stream1, rvec_gmm, el, bootMC, Kv_sim_h1, CVMv_sim_h1);
        QVrejvecs(j,:)       = QVrej_j;
        CVMrejvecs(j,:)      = CVMrej_j;
        QVrejvecs_bonf(j,:)  = QVbonf_j;
        CVMrejvecs_bonf(j,:) = CVMbonf_j;
    end
    results.(mod).gmm_ks       = mean(QVrejvecs,      1);
    results.(mod).gmm_cvm      = mean(CVMrejvecs,      1);
    results.(mod).gmm_ks_bonf  = mean(QVrejvecs_bonf, 1);
    results.(mod).gmm_cvm_bonf = mean(CVMrejvecs_bonf,1);
    fprintf('GMM done: %s\n', mod);
end
