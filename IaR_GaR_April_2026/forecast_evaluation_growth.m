%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FORECAST EVALUATION — GDP GROWTH  (multi-model comparison)
%%
%%  Authors: David Aikman, Rhys Bidder, Simon Lloyd, Giulia Mantoan,
%%           Simone Maso, Aditya Mori, Matthew Tong
%%
%%  Tests: KS · Rossi-Sekhposyan (2019) · Berkowitz (2001) ·
%%         Knueppel (2015) · Mitchell-Weale (2023) · Galvao-Mantoan-Mitchell
%%
%%  Compares: Fan Chart (BOE) · QR Skew-t (RASS) · BVAR
%%
%%  Notes:
%%    - Data are quarterly; h=1 = nowcast, h=2 = 1Q ahead, … h=13 = 12Q ahead.
%%      cfg.eval_horizons is in quarters (0..12); internal h_idx = k+1.
%%    - All three models are evaluated on exactly the same quarterly origins
%%      (eval_origins), the same actual g4rgdp realisations, and the same
%%      Covid exclusion window.
%%    - Actual data loaded ONCE from actual_gdp_yoy_OOS.mat; aligned to a
%%      single shared quarterly grid.
%%    - CRPS and left-tail wCRPS computed via deterministic quantile grid
%%      (Gneiting & Ranjan 2011, Eq. 15/17) — no Monte Carlo, no randomness.
%%      Left-tail weight: v(tau) = (1-tau)^2  [QW=5, GDP downside risk].
%%    - Fan chart uses TPN (Two-Piece Normal); RASS uses the 19-quantile grid
%%      via monotone interpolation; BVAR uses Normal(pred_mu, pred_sigma).
%%    - Covid origins (2020Q1–2022Q1) are set to NaN for all models.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
close all; clear; clc;
 
%% ── Paths ────────────────────────────────────────────────────────────────
scriptDir = fileparts(mfilename('fullpath'));
outDir    = fullfile(scriptDir, 'Outputs');
evalDir   = fullfile(outDir, 'forecast_evaluation');
if ~exist(evalDir,'dir'), mkdir(evalDir); end
 
%% ════════════════════════════════════════════════════════════════════════
%%  CONFIGURATION  ← edit here
%% ════════════════════════════════════════════════════════════════════════
 
% Horizons to evaluate (in quarters): 0 = nowcast, 1 = 1Q ahead, …, 12 = 3Y ahead
cfg.eval_horizons = 0:12;
 
% Common evaluation window  (pinned to fan chart availability)
cfg.eval_start = datetime('30-Sep-2007', 'InputFormat', 'dd-MMM-yyyy');
cfg.eval_end   = datetime('31-Mar-2024', 'InputFormat', 'dd-MMM-yyyy');
 
% Covid exclusion: origins in [covidStart, covidEnd] are set to NaN in all models.
cfg.covidStart = datetime(2020, 3, 31);   % Q1 2020
cfg.covidEnd   = datetime(2022, 3, 31);   % Q1 2022
 
% Fan chart settings
cfg.fanchart.start_date = datetime('30-Sep-2007', 'InputFormat', 'dd-MMM-yyyy');
cfg.fanchart.end_date   = datetime('31-Mar-2024', 'InputFormat', 'dd-MMM-yyyy');
cfg.fanchart.end_fcst   = datetime('30-Mar-2027', 'InputFormat', 'dd-MMM-yyyy');
cfg.fanchart.covid_date = datetime('30-Jun-2020', 'InputFormat', 'dd-MMM-yyyy');
 
% Root directory for boefsctdata_growth (needs Data\kcl_data\ subfolder)
cfg.paths.boe_root  = scriptDir;
cfg.paths.qr_pred_q = fullfile(outDir, 'predicted_quantiles_gdp_OOS.mat');
cfg.paths.qr_actual = fullfile(outDir, 'actual_gdp_yoy_OOS.mat');
cfg.paths.bvar      = fullfile(outDir, 'BVAR', 'BVAR_gdp_pred_q.mat');
 
%% ── Global settings ──────────────────────────────────────────────────────
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth', 1);
rng(0, 'twister');
 
addpath(fullfile(scriptDir, 'intermediate_codes'));
addpath(fullfile(scriptDir, 'functions'));
addpath(fullfile(scriptDir, 'functions', 'azzalini'));
addpath(fullfile(scriptDir, 'functions', 'CRPS'));
 
nHor    = numel(cfg.eval_horizons);
models  = {'fanchart', 'qr', 'bvar'};
nModels = numel(models);
 
% Shared deterministic quantile grid (999 points)
tau_grid = (1:999)' / 1000;
% Left-tail weight v(tau) = (1-tau)^2  [QW=5 — GDP downside / GaR focus]
v_left   = (1 - tau_grid).^2;
 
% CRPS and left-tail wCRPS helpers
crps_from_qgrid  = @(q, y) 2 * mean((double(y < q) - tau_grid) .* (q - y));
wcrps_from_qgrid = @(q, y) 2 * mean(v_left .* (double(y < q) - tau_grid) .* (q - y));
 
% Quantile levels used for sharpness interval widths
q_sharp_lev = [0.05, 0.10, 0.25, 0.35, 0.50, 0.65, 0.75, 0.90, 0.95];
nQL = numel(q_sharp_lev);
 
%% ════════════════════════════════════════════════════════════════════════
%%  SHARED ACTUAL DATA  (g4rgdp — loaded once, aligned to common grid)
%% ════════════════════════════════════════════════════════════════════════
fprintf('Building shared actual data (g4rgdp)...\n');
 
ad_shared = load(cfg.paths.qr_actual, 'actual_var', 'dateNumeric_full', 'idx_est');
gdp_series    = ad_shared.actual_var;           % full quarterly series
gdp_datenum   = ad_shared.dateNumeric_full;     % datenum for every row
gdp_ym        = floor(gdp_datenum / 30.4375);  % rough ym  — use datetime instead
% More precise: convert to datetime then compute ym
gdp_dt        = datetime(gdp_datenum, 'ConvertFrom', 'datenum');
gdp_ym        = year(gdp_dt)*12 + month(gdp_dt);
 
% Common quarterly evaluation grid
eval_origins = (cfg.eval_start : calmonths(3) : cfg.eval_end)';
nQOrig_eval  = numel(eval_origins);
eval_ym      = year(eval_origins)*12 + month(eval_origins);
 
% actualvar(ih, j) = g4rgdp realisation cfg.eval_horizons(ih) quarters
% ahead of eval_origins(j).   k quarters = 3k months in ym space.
actualvar = NaN(nHor, nQOrig_eval);
for j = 1:nQOrig_eval
    for ih = 1:nHor
        k      = cfg.eval_horizons(ih);
        tgt_ym = eval_ym(j) + 3*k;
        idx    = find(gdp_ym == tgt_ym, 1);
        if ~isempty(idx)
            actualvar(ih, j) = gdp_series(idx);
        end
    end
end
 
% Covid exclusion — two-dimensional:
%   (1) Origin-based: exclude any origin whose conditioning date is in the
%       Covid window (so no forecasts are made from Covid origins).
%   (2) Target-based: for each horizon h, also exclude any (origin, h) pair
%       whose TARGET date (origin + h quarters) falls in the Covid window.
%       Without this, pre-Covid origins with h>0 targets landing in the
%       Covid crash/rebound still score huge CRPS and inflate the results.
%
% Both masks extend 4 quarters beyond covidEnd to cover YoY base effects
% and align with the QR estimation window.
covidEnd_ext_eval = cfg.covidEnd + calquarters(4);   % Q1 2023
 
% (1) Origin mask
covid_mask = eval_origins >= cfg.covidStart & eval_origins <= covidEnd_ext_eval;
actualvar(:, covid_mask) = NaN;
 
% (2) Target mask: for each horizon, shift the Covid window backward by h
%     quarters so that (origin + h) ∈ Covid ⟹ NaN
for ih = 1:nHor
    k = cfg.eval_horizons(ih);
    % target date of origin j at horizon ih = eval_origins(j) + k quarters
    tgt_dates = eval_origins + calquarters(k);
    tgt_covid  = tgt_dates >= cfg.covidStart & tgt_dates <= covidEnd_ext_eval;
    actualvar(ih, tgt_covid) = NaN;
end
 
clearvars ad_shared gdp_series gdp_datenum gdp_dt;
 
%% ════════════════════════════════════════════════════════════════════════
%%  LOAD ALL MODELS AND COMPUTE PITs, CRPS, wCRPS, SHARPNESS
%%
%%  zinf_all{m}  — nHor × nQOrig_eval  PIT values
%%  crps_all{m}  — nHor × nQOrig_eval  uniform CRPS
%%  wcrps_all{m} — nHor × nQOrig_eval  left-tail wCRPS  [v=(1-tau)^2]
%%  qntl_all{m}  — nHor × nQOrig_eval × nQL  quantiles (sharpness)
%% ════════════════════════════════════════════════════════════════════════
 
zinf_all     = cell(nModels, 1);
crps_all     = cell(nModels, 1);
wcrps_all    = cell(nModels, 1);
qntl_all     = cell(nModels, 1);
originDT_all = cell(nModels, 1);
 
%% ── (1) Fan chart ────────────────────────────────────────────────────────
fprintf('Loading BOE fan chart...\n');
 
% boefsctdata_growth clears most of the workspace via clearvars -except.
% Workaround: save/restore workspace around the call.
start_date   = cfg.fanchart.start_date;
end_date     = cfg.fanchart.end_date;
end_fcst     = cfg.fanchart.end_fcst;
covid_date   = cfg.fanchart.covid_date;
fullFileName = cfg.paths.qr_actual;   % only needed by actualdata (not used below)
varnames     = {'g4rgdp'};
ctrynames    = {'UK'};
momentlist   = {'fcstmean','fcststdev','fcstskew'};
outputFolder = evalDir;
 
% boefsctdata_growth reads the Excel via [cd '\Data\kcl_data\...'],
% so we must be in scriptDir before calling it.
% After clearvars, only the exception variables survive (including outputFolder).
% We use outputFolder (= evalDir) to locate the backup and restore the workspace.
%cd(cfg.paths.boe_root);                          % set cwd BEFORE clearvars wipes cfg
save(fullfile(evalDir, 'ws_tmp_growth.mat'));     % save full workspace
boefsctdata_growth;   % produces mtestdata; clears all but exceptions
% scriptDir / evalDir / cfg are now gone — outputFolder survived
mtestdata_fc     = mtestdata;
load(fullfile(outputFolder, 'ws_tmp_growth.mat')); % restore full workspace
delete(fullfile(evalDir, 'ws_tmp_growth.mat'));
mtestdata = mtestdata_fc;
clear mtestdata_fc;
 
% Extract parameters: rows = horizons 1..13 (h=1 nowcast), cols = origins
modevar_fc = mtestdata(:, 1:3:end);   % 13 × nQOrig_fc
meanvar_fc = mtestdata(:, 2:3:end);
dispvar_fc = mtestdata(:, 3:3:end);
nQOrig_fc  = size(meanvar_fc, 2);
 
% Fan chart origin datetimes (quarterly from start_date)
fc_origins = (cfg.fanchart.start_date : calmonths(3) : ...
              cfg.fanchart.start_date + calmonths(3*(nQOrig_fc-1)))';
fc_ym      = year(fc_origins)*12 + month(fc_origins);
 
zinf_fc  = NaN(nHor, nQOrig_eval);
crps_fc  = NaN(nHor, nQOrig_eval);
wcrps_fc = NaN(nHor, nQOrig_eval);
qntl_fc  = NaN(nHor, nQOrig_eval, nQL);
 
for j = 1:nQOrig_eval
    fc_col = find(fc_ym == eval_ym(j), 1);
    if isempty(fc_col), continue; end
 
    for ih = 1:nHor
        % Fan chart row index: eval horizon k → h=k+1 (1-based)
        row = cfg.eval_horizons(ih) + 1;
        if row > size(meanvar_fc, 1), continue; end
        mn = meanvar_fc(row, fc_col);
        mo = modevar_fc(row, fc_col);
        dv = dispvar_fc(row, fc_col);
        if isnan(mn) || isnan(mo) || isnan(dv) || dv <= 0, continue; end
 
        if mn ~= mo
            % Asymmetric TPN
            gam = stdtogam(mn, mo, dv);
            sig = mom2g(dv, gam);
            [m_par, s1, s2] = momtopar(mn, mo, sig);
            p_left    = s1 / (s1 + s2);
            lm        = tau_grid <= p_left;
            q_grid    = NaN(size(tau_grid));
            q_grid( lm) = m_par + s1 .* norminv(tau_grid( lm) .* (s1+s2) ./ (2*s1));
            q_grid(~lm) = m_par - s2 .* norminv((1-tau_grid(~lm)) .* (s1+s2) ./ (2*s2));
            lm_s = q_sharp_lev <= p_left;
            q_s  = NaN(1, nQL);
            q_s( lm_s) = m_par + s1 .* norminv(q_sharp_lev( lm_s) .* (s1+s2) ./ (2*s1));
            q_s(~lm_s) = m_par - s2 .* norminv((1-q_sharp_lev(~lm_s)) .* (s1+s2) ./ (2*s2));
        else
            % Symmetric (normal) case
            q_grid = mn + dv .* norminv(tau_grid);
            q_s    = mn + dv .* norminv(q_sharp_lev);
            m_par  = mn;  s1 = dv;  s2 = dv;   % for ftp PIT
        end
 
        act = actualvar(ih, j);
        if isnan(act), continue; end
 
        % Sharpness set only for non-Covid (origin,horizon) pairs
        qntl_fc(ih, j, :) = q_s;
 
        zinf_fc(ih, j) = integral(@(y) ftp(y, m_par, s1, s2), -1e4, act, ...
                                   'RelTol', 1e-6, 'AbsTol', 1e-10);
        crps_fc(ih,j)  = crps_from_qgrid(q_grid, act);
        wcrps_fc(ih,j) = wcrps_from_qgrid(q_grid, act);
    end
end
 
zinf_all{1}     = zinf_fc;
crps_all{1}     = crps_fc;
wcrps_all{1}    = wcrps_fc;
qntl_all{1}     = qntl_fc;
originDT_all{1} = eval_origins;
 
clearvars modevar_fc meanvar_fc dispvar_fc fc_origins fc_ym fc_col ...
          zinf_fc crps_fc wcrps_fc qntl_fc mtestdata ...
          start_date end_date end_fcst covid_date fullFileName varnames ctrynames momentlist outputFolder ...
          gam sig ih j row m_par s1 s2 mn mo dv act p_left lm lm_s q_s q_grid nQOrig_fc;
 
%% ── (2) RASS Quantile Regression ─────────────────────────────────────────
%  The GDP RASS saves pred_q (nOrigins × 19 × 13) — 19 quantile levels per
%  origin and horizon.  PIT and CRPS are computed via monotone interpolation
%  of these 19 quantiles onto the 999-point tau_grid.  No skew-t fitting needed.
fprintf('Loading RASS model...\n');
 
qd = load(cfg.paths.qr_pred_q, 'pred_q');
ad = load(cfg.paths.qr_actual, 'actual_var', 'idx_est', 'dateNumeric_full');
 
pred_q_qr  = qd.pred_q;                % nOrig_qr × 19 × 13
quant_qr   = (0.05:0.05:0.95)';        % 19 × 1  (quantile levels in pred_q)
nOrig_qr   = size(pred_q_qr, 1);
 
qr_dt      = datetime(ad.dateNumeric_full(ad.idx_est : ad.idx_est + nOrig_qr - 1), ...
                       'ConvertFrom', 'datenum');
qr_ym      = year(qr_dt)*12 + month(qr_dt);
 
zinf_qr  = NaN(nHor, nQOrig_eval);
crps_qr  = NaN(nHor, nQOrig_eval);
wcrps_qr = NaN(nHor, nQOrig_eval);
qntl_qr  = NaN(nHor, nQOrig_eval, nQL);
 
for j = 1:nQOrig_eval
    t_idx = find(qr_ym == eval_ym(j), 1);
    if isempty(t_idx), continue; end
 
    for ih = 1:nHor
        k     = cfg.eval_horizons(ih);
        h_idx = k + 1;                           % pred_q h-index: 1=nowcast
        if h_idx > size(pred_q_qr, 3), continue; end
 
        q_vals = squeeze(pred_q_qr(t_idx, :, h_idx))';   % 19 × 1
        if all(isnan(q_vals)), continue; end
 
        act = actualvar(ih, j);
        if isnan(act), continue; end
 
        % Sharpness set only for non-Covid (origin,horizon) pairs
        qntl_qr(ih, j, :) = interp1(quant_qr, q_vals, q_sharp_lev, ...
                                     'linear', 'extrap');
 
        % PIT = F(act) = tau such that Q(tau) = act
        %       interp1 on (q_vals → quant_qr) inverts the quantile function
        pit = interp1(q_vals, quant_qr, act, 'linear', 'extrap');
        zinf_qr(ih, j) = min(max(pit, 0), 1);
 
        % Expand 19-quantile grid to 999 points for CRPS
        q_fine = interp1(quant_qr, q_vals, tau_grid, 'linear', 'extrap');
 
        crps_qr(ih, j)  = crps_from_qgrid(q_fine, act);
        wcrps_qr(ih, j) = wcrps_from_qgrid(q_fine, act);
    end
end
 
zinf_all{2}     = zinf_qr;
crps_all{2}     = crps_qr;
wcrps_all{2}    = wcrps_qr;
qntl_all{2}     = qntl_qr;
originDT_all{2} = eval_origins;
 
clearvars qd ad pred_q_qr quant_qr nOrig_qr qr_dt qr_ym ...
          zinf_qr crps_qr wcrps_qr qntl_qr ...
          ih j k h_idx t_idx q_vals act pit q_fine;
 
%% ── (3) BVAR  (Normal predictive — mean and std from Gibbs draws) ────────
fprintf('Loading BVAR...\n');
 
bd = load(cfg.paths.bvar, 'pred_mu', 'pred_sigma', 'dateNumeric_est');
pred_mu_bv    = bd.pred_mu;      % nOrig_bv × 13
pred_sigma_bv = bd.pred_sigma;   % nOrig_bv × 13
 
bv_dt = datetime(bd.dateNumeric_est, 'ConvertFrom', 'datenum')';
bv_ym = year(bv_dt)*12 + month(bv_dt);
 
zinf_bv  = NaN(nHor, nQOrig_eval);
crps_bv  = NaN(nHor, nQOrig_eval);
wcrps_bv = NaN(nHor, nQOrig_eval);
qntl_bv  = NaN(nHor, nQOrig_eval, nQL);
 
for j = 1:nQOrig_eval
    t_idx = find(bv_ym == eval_ym(j), 1);
    if isempty(t_idx), continue; end
 
    for ih = 1:nHor
        k     = cfg.eval_horizons(ih);
        h_idx = k + 1;
        if h_idx > size(pred_mu_bv, 2), continue; end
 
        mu_bv    = pred_mu_bv(t_idx, h_idx);
        sigma_bv = pred_sigma_bv(t_idx, h_idx);
        if isnan(mu_bv) || isnan(sigma_bv) || sigma_bv <= 0, continue; end
 
        act = actualvar(ih, j);
        if isnan(act), continue; end
 
        % Sharpness set only for non-Covid (origin,horizon) pairs
        qntl_bv(ih, j, :) = norminv(q_sharp_lev, mu_bv, sigma_bv);
 
        zinf_bv(ih, j) = normcdf(act, mu_bv, sigma_bv);
 
        q_grid = norminv(tau_grid, mu_bv, sigma_bv);
        crps_bv(ih, j)  = crps_from_qgrid(q_grid, act);
        wcrps_bv(ih, j) = wcrps_from_qgrid(q_grid, act);
    end
end
 
zinf_all{3}     = zinf_bv;
crps_all{3}     = crps_bv;
wcrps_all{3}    = wcrps_bv;
qntl_all{3}     = qntl_bv;
originDT_all{3} = eval_origins;
 
clearvars bd pred_mu_bv pred_sigma_bv bv_dt bv_ym ...
          zinf_bv crps_bv wcrps_bv qntl_bv ...
          ih j k h_idx t_idx mu_bv sigma_bv act q_grid;
 
 
%% ════════════════════════════════════════════════════════════════════════
%%  RUN TESTS FOR EACH MODEL
%% ════════════════════════════════════════════════════════════════════════
 
results = struct();
for m = 1:nModels
    results.(models{m}).ksinf        = NaN(nHor, 1);
    results.(models{m}).ksinfpv      = NaN(nHor, 1);
    results.(models{m}).rs_ks_test   = NaN(nHor, 1);
    results.(models{m}).rs_cvm_test  = NaN(nHor, 1);
    results.(models{m}).rs_ks_logic  = NaN(nHor, 3);
    results.(models{m}).rs_cvm_logic = NaN(nHor, 3);
    results.(models{m}).bert1        = NaN(nHor, 1);
    results.(models{m}).bert2        = NaN(nHor, 1);
    results.(models{m}).berK1        = NaN(nHor, 1);
    results.(models{m}).berK2        = NaN(nHor, 1);
    results.(models{m}).knueppel_stat= NaN(nHor, 1);
    results.(models{m}).knueppel_pval= NaN(nHor, 1);
    results.(models{m}).MW_stat      = NaN(nHor, 1);
    results.(models{m}).MW_pval      = NaN(nHor, 1);
    results.(models{m}).crps_mean    = NaN(nHor, 1);
    results.(models{m}).wcrps_mean   = NaN(nHor, 1);   % left-tail wCRPS [v=(1-tau)^2]
end
 
lags     = -1;
prewhite =  0;
z_u = 0.95;  z_l = 0.05;  df = 5;
rvec_rs = linspace(0, 1, 1000);
 
for m = 1:nModels
    mod  = models{m};
    zinf = zinf_all{m};
    fprintf('\nRunning tests: %s ...\n', mod);
 
    KS_vec  = zeros(nHor, 1);
    CVM_vec = zeros(nHor, 1);
 
    for i = 1:nHor
        z_row = zinf(i, :)';
        z_row = z_row(~isnan(z_row));
        if numel(z_row) < 5, continue; end
 
        % KS
        [results.(mod).ksinf(i), results.(mod).ksinfpv(i)] = kstestu(z_row);
 
        % Rossi-Sekhposyan (2019)
        [~, stat_i, logic_i, KS_vec(i), CVM_vec(i)] = rs_test(z_row, rvec_rs);
        results.(mod).rs_ks_test(i)    = stat_i.ks_stat;
        results.(mod).rs_cvm_test(i)   = stat_i.cvm_stat;
        results.(mod).rs_ks_logic(i,:)  = logic_i.ks_array;
        results.(mod).rs_cvm_logic(i,:) = logic_i.cvm_array;
 
        % Berkowitz (2001)
        [results.(mod).bert1(i), results.(mod).bert2(i), ...
         results.(mod).berK1(i), results.(mod).berK2(i)] = berk(zinf(i, :));
 
        % Knueppel (2015)
        [results.(mod).knueppel_stat(i), results.(mod).knueppel_pval(i)] = ...
            alpha0_1234_NW(zinf(i,:), lags, prewhite);
 
        % Mitchell-Weale (2023)
        stat_u = MW_alpha0_1234_NW(zinf(i,:), lags, prewhite, z_l, z_u);
        stat_f = freqTestInCensoredRegion(zinf(i,:), lags, prewhite, z_l, z_u);
        results.(mod).MW_stat(i) = stat_u + stat_f;
        results.(mod).MW_pval(i) = 1 - chi2cdf(results.(mod).MW_stat(i), df);
 
        % CRPS mean (uniform)
        crps_row = crps_all{m}(i, :);
        crps_row = crps_row(~isnan(crps_row));
        if ~isempty(crps_row)
            results.(mod).crps_mean(i) = mean(crps_row);
        end
 
        % wCRPS mean (left-tail weighted, v=(1-tau)^2)
        wcrps_row = wcrps_all{m}(i, :);
        wcrps_row = wcrps_row(~isnan(wcrps_row));
        if ~isempty(wcrps_row)
            results.(mod).wcrps_mean(i) = mean(wcrps_row);
        end
    end
 
    results.(mod).KS_vec  = KS_vec;
    results.(mod).CVM_vec = CVM_vec;
end
 
%% Galvao-Mantoan-Mitchell  (computationally intensive — runs per model)
MC     = 1000;
bootMC = 1000;
rng(bootMC, 'twister');
rvec_gmm = 0:0.001:1;

for m = 1:nModels
    mod  = models{m};
    zinf = zinf_all{m};
    z    = zinf';
    P    = size(z, 1);
    el   = floor(P^(1/4));
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
            size_statistic_h2(z, KS, CVM, Hz, stream1, rvec_gmm, el, bootMC);
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
 
%% ════════════════════════════════════════════════════════════════════════
%%  SAVE COMPARISON RESULTS
%% ════════════════════════════════════════════════════════════════════════
 
horiz    = cfg.eval_horizons(:);
filename = fullfile(evalDir, 'comparison_fcst_eval_growth.xlsx');
 
buildTab = @(field) array2table( ...
    [horiz, cell2mat(cellfun(@(m) results.(m).(field), models(:)', 'UniformOutput', false))], ...
    'VariableNames', [{'horizon'}, models(:)']);
 
writetable(buildTab('ksinf'),         filename, 'Sheet', 'KS_stat');
writetable(buildTab('ksinfpv'),       filename, 'Sheet', 'KS_pval');
writetable(buildTab('rs_ks_test'),    filename, 'Sheet', 'RS_KS_stat');
writetable(buildTab('rs_cvm_test'),   filename, 'Sheet', 'RS_CVM_stat');
writetable(buildTab('bert1'),         filename, 'Sheet', 'Berk_rho0_stat');
writetable(buildTab('bert2'),         filename, 'Sheet', 'Berk_rhohat_stat');
writetable(buildTab('berK1'),         filename, 'Sheet', 'Berk_rho0_pval');
writetable(buildTab('berK2'),         filename, 'Sheet', 'Berk_rhohat_pval');
writetable(buildTab('knueppel_stat'), filename, 'Sheet', 'Knueppel_stat');
writetable(buildTab('knueppel_pval'), filename, 'Sheet', 'Knueppel_pval');
writetable(buildTab('MW_stat'),       filename, 'Sheet', 'MW_stat');
writetable(buildTab('MW_pval'),       filename, 'Sheet', 'MW_pval');
writetable(buildTab('crps_mean'),     filename, 'Sheet', 'CRPS_mean');
writetable(buildTab('wcrps_mean'),    filename, 'Sheet', 'wCRPS_left_mean');
 
gmm_ks_mat  = cell2mat(cellfun(@(m) results.(m).gmm_ks(:)',  models(:)', 'UniformOutput', false)')';
gmm_cvm_mat = cell2mat(cellfun(@(m) results.(m).gmm_cvm(:)', models(:)', 'UniformOutput', false)')';
gmm_tab = table(models', gmm_ks_mat, gmm_cvm_mat, ...
    'VariableNames', {'model','GMM_KS_std_bonf_w_wi','GMM_CVM_std_bonf_w_wi'});
writetable(gmm_tab, filename, 'Sheet', 'GMM_summary');
 
%% ════════════════════════════════════════════════════════════════════════
%%  PRINT SUMMARY TABLES  (one per model)
%% ════════════════════════════════════════════════════════════════════════
 
model_labels = {'Fan Chart (BOE)', 'QR Skew-t', 'BVAR'};
 
sel_h   = [0 1 2 3 4 8 12];
sel_idx = arrayfun(@(h) find(cfg.eval_horizons == h, 1), sel_h);
% sel_h and sel_idx used both in the plain text printer and the LaTeX printer
 
frcode = @(v) char(70*(~isnan(v) & v < 0.5) + ...
                   82*(~isnan(v) & v >= 0.5) + ...
                   45*isnan(v));
 
hdr = ['  h   | Knu pval  |  MW pval  | CRPS mean | wCRPS(L)' ...
       ' | CvM 10% | CvM  5% | CvM  1%'];
sep = repmat('-', 1, numel(hdr));
 
for m = 1:nModels
    mod = models{m};
    fprintf('\n%s\n', sep);
    fprintf('  %s\n', model_labels{m});
    fprintf('%s\n', sep);
    fprintf('%s\n', hdr);
    fprintf('%s\n', sep);
    for ii = 1:numel(sel_idx)
        i   = sel_idx(ii);
        h   = cfg.eval_horizons(i);
        cvm = results.(mod).rs_cvm_logic(i, :);
        fprintf('  %2d  |  %7.4f  |  %7.4f  | %9.4f | %8.4f |   %s     |   %s     |   %s\n', ...
            h, ...
            results.(mod).knueppel_pval(i), ...
            results.(mod).MW_pval(i), ...
            results.(mod).crps_mean(i), ...
            results.(mod).wcrps_mean(i), ...
            frcode(cvm(3)), frcode(cvm(2)), frcode(cvm(1)));
    end
    fprintf('%s\n', sep);
end
 
%% ════════════════════════════════════════════════════════════════════════
%%  LATEX TABLE  — p-values + CRPS + wCRPS, all models, selected horizons
%%
%%  Paste directly into Overleaf.  Significance stars on p-values:
%%    *** p < 0.01  |  ** p < 0.05  |  * p < 0.10
%% ════════════════════════════════════════════════════════════════════════
 
pstar = @(p) [repmat('*', 1, (p<0.10) + (p<0.05) + (p<0.01))];
 
% Column header: one column per model
col_spec = ['l', repmat('c', 1, nModels)];
col_hdrs = strjoin(cellfun(@(l) ['\textbf{', strrep(l,' ','\;'), '}'], ...
                   model_labels, 'UniformOutput', false), ' & ');
 
latex_file = fullfile(evalDir, 'table_fcst_eval_growth.tex');
fid = fopen(latex_file, 'w');
 
fprintf(fid, '%% GDP forecast evaluation table — auto-generated\n');
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{GDP Growth Forecast Evaluation}\n');
fprintf(fid, '\\label{tab:fcst_eval_growth}\n');
fprintf(fid, '\\small\n');
fprintf(fid, '\\begin{tabular}{%s}\n', col_spec);
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Horizon & %s \\\\\n', col_hdrs);
fprintf(fid, '\\midrule\n');
 
% Block printer (inline to avoid nested-function scoping issues in MATLAB scripts)
tex_blocks = { ...
    'Knueppel (2015) p-value', 'knueppel_pval',  true;  ...
    'Mitchell-Weale p-value',  'MW_pval',        true;  ...
    'CRPS (mean)',             'crps_mean',      false; ...
    'wCRPS left tail (mean)',  'wcrps_mean',     false  ...
};
 
for tb = 1:size(tex_blocks, 1)
    blk_label  = tex_blocks{tb, 1};
    blk_field  = tex_blocks{tb, 2};
    blk_ispval = tex_blocks{tb, 3};
    if tb > 1, fprintf(fid, '\\addlinespace\n'); end
    fprintf(fid, '\\multicolumn{%d}{l}{\\textit{%s}} \\\\\n', nModels+1, blk_label);
    for ii = 1:numel(sel_idx)
        si   = sel_idx(ii);
        h_ii = sel_h(ii);
        row  = sprintf('$h=%d$', h_ii);
        for m = 1:nModels
            v = results.(models{m}).(blk_field)(si);
            if isnan(v)
                row = [row, ' & ---'];
            elseif blk_ispval
                row = [row, sprintf(' & %.3f%s', v, pstar(v))];
            else
                row = [row, sprintf(' & %.3f', v)];
            end
        end
        fprintf(fid, '%s \\\\\n', row);
    end
end
 
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\begin{tablenotes}\\small\n');
fprintf(fid, '\\item Notes: $h$ is horizon in quarters (0 = nowcast). ');
fprintf(fid, 'p-values: $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$. ');
fprintf(fid, 'wCRPS uses left-tail weight $v(\\tau)=(1-\\tau)^2$. ');
fprintf(fid, 'Covid origins (2020Q1--2022Q1) excluded.\n');
fprintf(fid, '\\end{tablenotes}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
 
% Also print to console
fprintf('\n\n%% ── LaTeX table ─────────────────────────────────────────\n');
fid2 = fopen(latex_file, 'r');
while ~feof(fid2), fprintf('%s\n', fgetl(fid2)); end
fclose(fid2);
fprintf('\nLaTeX table saved to %s\n', latex_file);
 
%% ════════════════════════════════════════════════════════════════════════
%%  SHARPNESS DIAGNOSTIC
%%
%%  Quantile indices in qntl_all{m}(ih, j, :):
%%    1=q05  2=q10  3=q25  4=q35  5=q50  6=q65  7=q75  8=q90  9=q95
%% ════════════════════════════════════════════════════════════════════════
fprintf('\nComputing sharpness widths...\n');
 
h_vec_sharp  = cfg.eval_horizons(:);
% 2×2 layout: Central 30% | Central 50% | Bottom 20% | Top 20%
sharp_labels = {'Central 30%','Central 50%','Bottom 20%','Top 20%'};
nW = numel(sharp_labels);
 
mu_w = cell(nModels, 1);
sd_w = cell(nModels, 1);
 
for m = 1:nModels
    Q  = qntl_all{m};   % nHor × nQOrig_eval × 9
    % Quantile indices: 1=q05 2=q10 3=q25 4=q35 5=q50 6=q65 7=q75 8=q90 9=q95
    mu_w{m} = NaN(nHor, nW);
    sd_w{m} = NaN(nHor, nW);
    for ih = 1:nHor
        widths = [ squeeze(Q(ih,:,6)) - squeeze(Q(ih,:,4));   % w30
                   squeeze(Q(ih,:,7)) - squeeze(Q(ih,:,3));   % w50
                   squeeze(Q(ih,:,3)) - squeeze(Q(ih,:,1));   % wb20
                   squeeze(Q(ih,:,9)) - squeeze(Q(ih,:,7)) ]; % wt20
        mu_w{m}(ih,:) = median(widths, 2, 'omitnan')';
        sd_w{m}(ih,:) = std(widths, 0, 2, 'omitnan')';
    end
end
 
% Print sharpness table
fprintf('\n%s\n', repmat('─',1,70));
fprintf('  SHARPNESS — Median interval widths by horizon\n');
fprintf('%s\n', repmat('─',1,70));
fprintf('  h  | %s\n', strjoin(cellfun(@(s) sprintf('%11s',s), sharp_labels, 'UniformOutput',false),' | '));
for m = 1:nModels
    fprintf('\n  %s\n', model_labels{m});
    fprintf('  %s\n', repmat('-',1,68));
    for ii = 1:numel(sel_idx)
        i = sel_idx(ii);
        h = cfg.eval_horizons(i);
        fprintf('  %2d |', h);
        for iw = 1:nW
            fprintf('  %9.3f |', mu_w{m}(i, iw));
        end
        fprintf('\n');
    end
end
fprintf('%s\n', repmat('─',1,70));
 
% ── Shared style ─────────────────────────────────────────────────────────
colors_sh  = [0.12 0.47 0.71; 0.84 0.15 0.16; 0.17 0.63 0.17];
fills_sh   = [0.70 0.80 1.00; 1.00 0.80 0.80; 0.80 1.00 0.80];
lssh       = {'-', '--', ':'};
markers_sh = {'o', 's', '^'};
lwsh       = 1.8;
h_vec      = cfg.eval_horizons(:);
 
%% ── CHART 1: Sharpness — 2×2 ─────────────────────────────────────────────
fig_sh = figure('Name','Sharpness — GDP Growth','NumberTitle','off', ...
                'Color','w','Position',[120 120 900 750]);
tl_sh  = tiledlayout(fig_sh, 2, 2, 'TileSpacing','compact','Padding','loose');
sgtitle(fig_sh, 'GDP Growth — Forecast Sharpness (median width ± 1 SD)', ...
        'FontSize', 12, 'FontWeight','bold');
 
legH_sh = gobjects(nModels, 1);
for iw = 1:nW
    ax = nexttile(tl_sh, iw);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for m = 1:nModels
        mu_v = mu_w{m}(:, iw);
        sd_v = sd_w{m}(:, iw);
        xx   = [h_vec; flipud(h_vec)];
        yy   = [mu_v - sd_v; flipud(mu_v + sd_v)];
        patch('Parent',ax,'XData',xx,'YData',yy, ...
              'FaceColor',fills_sh(m,:),'EdgeColor','none', ...
              'FaceAlpha',0.40,'HandleVisibility','off');
        hh = plot(ax, h_vec, mu_v, lssh{m}, ...
                  'Color',colors_sh(m,:),'LineWidth',lwsh, ...
                  'Marker',markers_sh{m},'MarkerSize',4, ...
                  'MarkerFaceColor',colors_sh(m,:), ...
                  'DisplayName',model_labels{m});
        if iw == 1, legH_sh(m) = hh; end
    end
    xlabel(ax,'Horizon (quarters ahead)','FontSize',9);
    ylabel(ax,'Width (pp)','FontSize',9);
    title(ax, sharp_labels{iw},'FontSize',10,'FontWeight','bold');
    set(ax,'XTick',h_vec,'Box','on','FontSize',9);
    hold(ax,'off');
end
lgd_sh = legend(legH_sh, model_labels, 'Orientation','horizontal','NumColumns',3);
try, lgd_sh.Layout.Tile = 'south'; catch
    set(lgd_sh,'Position',[0.25 0.01 0.50 0.04]); end
exportgraphics(fig_sh, fullfile(evalDir,'sharpness_fcst_eval_growth.png'), 'Resolution',300);
fprintf('\nSharpness chart saved.\n');
 
%% ── CHART 2: Calibration p-values — 1×2 (Knueppel + Mitchell-Weale) ─────
fig_pv = figure('Name','Calibration Tests — GDP Growth','NumberTitle','off', ...
                'Color','w','Position',[80 80 900 420]);
tl_pv  = tiledlayout(fig_pv, 1, 2, 'TileSpacing','compact','Padding','loose');
sgtitle(fig_pv,'GDP Growth — Calibration Tests (p-values)', ...
        'FontSize',12,'FontWeight','bold');
 
pval_fields = {'knueppel_pval','MW_pval'};
pval_titles = {'Knueppel (2015) p-value','Mitchell-Weale (2023) p-value'};
legH_pv = gobjects(nModels,1);
 
for t = 1:2
    ax = nexttile(tl_pv, t);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for m = 1:nModels
        hh = plot(ax, h_vec, results.(models{m}).(pval_fields{t}), ...
             'Color',colors_sh(m,:),'LineStyle',lssh{m},'LineWidth',lwsh, ...
             'Marker',markers_sh{m},'MarkerSize',4,'MarkerFaceColor',colors_sh(m,:), ...
             'DisplayName',model_labels{m});
        if t == 1, legH_pv(m) = hh; end
    end
    yline(ax,0.10,'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',1.0, ...
          'Label','10%','LabelVerticalAlignment','bottom','HandleVisibility','off');
    yline(ax,0.05,'Color',[0.2 0.2 0.2],'LineStyle','--','LineWidth',1.0, ...
          'Label','5%','LabelVerticalAlignment','bottom','HandleVisibility','off');
    ylim(ax,[0 1]);
    xlabel(ax,'Horizon (quarters ahead)','FontSize',9);
    title(ax, pval_titles{t},'FontSize',10,'FontWeight','bold');
    set(ax,'XTick',h_vec,'Box','on','FontSize',9);
    hold(ax,'off');
end
lgd_pv = legend(legH_pv, model_labels,'Orientation','horizontal','NumColumns',3);
try, lgd_pv.Layout.Tile = 'south'; catch
    set(lgd_pv,'Position',[0.20 0.01 0.60 0.04]); end
exportgraphics(fig_pv, fullfile(evalDir,'pvalues_fcst_eval_growth.png'), 'Resolution',300);
fprintf('\nCalibration p-value chart saved.\n');
 
%% ── CHART 3: Cumulative CRPS & wCRPS — 2×2 ──────────────────────────────
cum_h       = [1, 4];
cum_h_names = {'1Q ahead (h=1)', '4Q ahead (h=4)'};
 
fig_sc = figure('Name','Cumulative Scores — GDP Growth','NumberTitle','off', ...
                'Color','w','Position',[100 100 1000 750]);
tl_sc  = tiledlayout(fig_sc, 2, 2, 'TileSpacing','compact','Padding','loose');
sgtitle(fig_sc,'GDP Growth — Cumulative CRPS & wCRPS (left tail)', ...
        'FontSize',12,'FontWeight','bold');
 
legH_sc = gobjects(nModels,1);
 
% Row 1: Cumulative CRPS at h=1 and h=4
for ip = 1:2
    ax_c = nexttile(tl_sc, ip);
    hold(ax_c,'on'); box(ax_c,'on'); grid(ax_c,'on');
    ih_row = find(cfg.eval_horizons == cum_h(ip), 1);
    for m = 1:nModels
        score_row = crps_all{m}(ih_row, :);
        valid = ~isnan(score_row);
        hh = plot(ax_c, eval_origins, cumsum(score_row .* valid), ...
             'Color',colors_sh(m,:),'LineStyle',lssh{m},'LineWidth',lwsh, ...
             'DisplayName',model_labels{m});
        if ip == 1, legH_sc(m) = hh; end
    end
    xlabel(ax_c,'Evaluation origin','FontSize',9);
    ylabel(ax_c,'Cumulative CRPS','FontSize',9);
    title(ax_c, sprintf('Cum. CRPS — %s', cum_h_names{ip}), ...
          'FontSize',10,'FontWeight','bold');
    set(ax_c,'FontSize',9);
    datetick(ax_c,'x','yyyy','keepticks','keeplimits'); xtickangle(ax_c,45);
    hold(ax_c,'off');
end
 
% Row 2: Cumulative wCRPS (left tail) at h=1 and h=4
for ip = 1:2
    ax_w = nexttile(tl_sc, ip+2);
    hold(ax_w,'on'); box(ax_w,'on'); grid(ax_w,'on');
    ih_row = find(cfg.eval_horizons == cum_h(ip), 1);
    for m = 1:nModels
        score_row = wcrps_all{m}(ih_row, :);
        valid = ~isnan(score_row);
        plot(ax_w, eval_origins, cumsum(score_row .* valid), ...
             'Color',colors_sh(m,:),'LineStyle',lssh{m},'LineWidth',lwsh, ...
             'DisplayName',model_labels{m});
    end
    xlabel(ax_w,'Evaluation origin','FontSize',9);
    ylabel(ax_w,'Cumulative wCRPS','FontSize',9);
    title(ax_w, sprintf('Cum. wCRPS (left tail) — %s', cum_h_names{ip}), ...
          'FontSize',10,'FontWeight','bold');
    set(ax_w,'FontSize',9);
    datetick(ax_w,'x','yyyy','keepticks','keeplimits'); xtickangle(ax_w,45);
    hold(ax_w,'off');
end
 
lgd_sc = legend(legH_sc, model_labels,'Orientation','horizontal','NumColumns',3);
try, lgd_sc.Layout.Tile = 'south'; catch
    set(lgd_sc,'Position',[0.20 0.01 0.60 0.04]); end
exportgraphics(fig_sc, fullfile(evalDir,'scores_fcst_eval_growth.png'), 'Resolution',300);
fprintf('\nCumulative scores chart saved.\n');
 
fprintf('\n── DONE: comparison saved to %s ──\n', filename);