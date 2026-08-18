function [result_sup, result_bonf, result_h] = CVfinalbootstrsuptest_shifted_h1rankfix(el, bootMC, pit, rvec, w, wi, stream1, Kv_sim_h1, CVMv_sim_h1)
% 2026-08-13: real-data counterpart of CVfinalbootstrsuptest_shifted2_h1rankfix.m.
% That version's cross-horizon correlation weighting (the "shifted2"
% eta/G construction) is derived FROM the simulated DGP's known theta
% parameter -- it hard-codes the IMA(1) correlation decay S(h) =
% sum_k theta^(2k), which has no meaningful value for real data with an
% unknown true dependence structure. This version instead builds on
% CVfinalbootstrsuptest_shifted.m's FULL-SHARING construction (same eta
% value at the index-aligned shift across horizons, correlation=1,
% no theta needed) -- a defensible, model-free choice that doesn't assume
% zero correlation (like the original, pre-shift bootstrap) or a specific
% wrong decay parameter (like shifted2 applied off-DGP).
%
% Horizon 1's marginal/Sup Test contribution is corrected the same way as
% CVfinalbootstrsuptest_shifted2_h1rankfix.m: KSv(:,1)/CVMv(:,1) are
% rank-preserving quantile-mapped onto the closed-form Brownian-bridge
% distribution (Kv_sim_h1, CVMv_sim_h1 -- precompute once via
% closedform_h1_cv_anchored.m, using the SAME rvec grid used here, and
% pass in), which is exactly R&S(2019)'s own treatment of h=1 and is
% fully generic -- it only relies on h=1 having no MA dependence under
% the null of correct specification, true for any model, not just the
% simulated DGP.

H=size(pit,2);
KSv = zeros(bootMC,H); CVMv = zeros(bootMC,H);

P = size(pit,1);
size_rvec = length(rvec);
n_windows = P - el + 1;

cs = cell(1,H);
for zvec = 1:H
    emp_cdf = (pit(:,zvec) <= rvec);
    mean_emp_cdf = mean(emp_cdf,1);
    centered = emp_cdf - mean_emp_cdf;
    cs{zvec} = [zeros(1,size_rvec); cumsum(centered,1)];
end

blocksums = cell(1,H);
for zvec = 1:H
    blocksums{zvec} = cs{zvec}(el+1:el+n_windows,:) - cs{zvec}(1:n_windows,:);
end

ETA_shared = (1/sqrt(el)) * randn(n_windows + H, bootMC);

for zvec = 1:H
    ETA_h = ETA_shared(zvec+1 : zvec+n_windows, :);   % n_windows x bootMC

    K_star_all = (P^(-1/2)) * (ETA_h' * blocksums{zvec});
    KSv(:,zvec) = max(abs(K_star_all),[],2);
    CVMv(:,zvec) = mean(K_star_all.^2,2);
end

% ---- rank-preserving quantile mapping for horizon 1 (see
% CVfinalbootstrsuptest_shifted2_h1rankfix.m header for the full
% rationale) ----
N_sim = numel(Kv_sim_h1);
[~, sortIdxK] = sort(KSv(:,1), 'ascend');
ranksK = zeros(bootMC,1); ranksK(sortIdxK) = 1:bootMC;
mapIdxK = max(1, min(N_sim, round(ranksK/bootMC * N_sim)));
KSv(:,1) = Kv_sim_h1(mapIdxK);

[~, sortIdxM] = sort(CVMv(:,1), 'ascend');
ranksM = zeros(bootMC,1); ranksM(sortIdxM) = 1:bootMC;
mapIdxM = max(1, min(N_sim, round(ranksM/bootMC * N_sim)));
CVMv(:,1) = CVMv_sim_h1(mapIdxM);

KSv_sup_max = max(KSv,[],2);
CVMv_sup_max = max(CVMv,[],2);

KSv_sup_w = max(KSv.*w,[],2);
CVMv_sup_w = max(CVMv.*w,[],2);

KSv_sup_wi = max(KSv.*wi,[],2);
CVMv_sup_wi = max(CVMv.*wi,[],2);

KSv = sort(KSv,'ascend');       cvKv_h = KSv(bootMC*0.95,:);
CVMv = sort(CVMv,'ascend');     cvMv_h = CVMv(bootMC*0.95,:);

KSv_sup_max = sort(KSv_sup_max,'ascend');       cvKv_max = KSv_sup_max(bootMC*0.95);
CVMv_sup_max = sort(CVMv_sup_max,'ascend');     cvMv_max = CVMv_sup_max(bootMC*0.95);

KSv_sup_w = sort(KSv_sup_w,'ascend');       cvKv_w = KSv_sup_w(bootMC*0.95);
CVMv_sup_w = sort(CVMv_sup_w,'ascend');     cvMv_w = CVMv_sup_w(bootMC*0.95);

KSv_sup_wi = sort(KSv_sup_wi,'ascend');       cvKv_wi = KSv_sup_wi(bootMC*0.95);
CVMv_sup_wi = sort(CVMv_sup_wi,'ascend');     cvMv_wi = CVMv_sup_wi(bootMC*0.95);

KSv_bonf = sort(KSv,'ascend');       cvKv_bonf = KSv_bonf(round(bootMC*(1-(0.05/H)),0),:);
CVMv_bonf = sort(CVMv,'ascend');     cvMv_bonf = CVMv_bonf(round(bootMC*(1-(0.05/H)),0),:);

result_sup = [cvKv_max cvMv_max; cvKv_w cvMv_w; cvKv_wi cvMv_wi];
result_bonf=[cvKv_bonf;cvMv_bonf];
result_h=[cvKv_h; cvMv_h];
end
