%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FORECAST EVALUATION — CPI INFLATION  (multi-model comparison)
%%
%%  Authors: David Aikman, Rhys Bidder, Simon Lloyd, Giulia Mantoan,
%%           Simone Maso, Aditya Mori, Matthew Tong
%%
%%  Tests: KS · Rossi-Sekhposyan (2019) · Berkowitz (2001) ·
%%         Knueppel (2015) · Mitchell-Weale (2023) · Galvao-Mantoan-Mitchell
%%
%%  Compares: Fan Chart (BOE) · QR Skew-t (RASS) · QR Semi-param (RASS) · BVAR
%%  Results are saved side-by-side for each test statistic.
%%
%%  Notes:
%%    - Evaluation is on QUARTERLY origins (months 3,6,9,12) and
%%      QUARTERLY horizons (0,3,6,...,36 months ahead = 0..12 quarters ahead).
%%    - No Covid exclusion applied to inflation.
%%    - Actual data: g4cpi (ONS headline CPI) loaded ONCE from the monthly
%%      Excel file, aligned to a single shared evaluation grid (eval_origins).
%%    - All three models are evaluated on exactly the same OOS dates.
%%    - g4Infl in GaRDataRaw.xlsx is CPIH (a different measure) and is NOT used.
%%    - CRPS and left-tail wCRPS computed via deterministic quantile grid
%%      (Gneiting & Ranjan 2011, Eq. 15/17) — no paretotails, no randomness.
%%      Left-tail weight: v(tau) = (1-tau)^2  [QW=5 in qwps notation].
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

% Horizons to evaluate: 0 = nowcast, 1 = 1Q ahead, …, 12 = 3Y ahead
cfg.eval_horizons = 0:12;

% No Covid exclusion for inflation
cfg.covid_exclude = false;

% Fan chart settings  (eval window is pinned to these dates)
cfg.fanchart.start_date = datetime('31-Mar-2010', 'InputFormat', 'dd-MMM-yyyy');
cfg.fanchart.end_date   = datetime('30-Sep-2022', 'InputFormat', 'dd-MMM-yyyy');
cfg.fanchart.end_fcst   = datetime('30-Sep-2025', 'InputFormat', 'dd-MMM-yyyy');

% Evaluation window is pinned to the fan chart origin window
% so all three models are always evaluated over exactly the same dates.
cfg.eval_start = cfg.fanchart.start_date;
cfg.eval_end   = cfg.fanchart.end_date;

% Root directory that boefsctdata uses via cd() to locate
% Data\cpi_infl_projection_parameters_mpc.xlsx
cfg.paths.boe_root    = scriptDir;
cfg.paths.monthly     = fullfile(scriptDir, 'IaRDataRaw_monthly_M.xlsx');
cfg.paths.qr_sktparam = fullfile(outDir, 'sktparam', 'skewtparam_inflation_OOS.mat');
cfg.paths.qr_actual   = fullfile(outDir, 'actual_inflation_mom_OOS.mat');
cfg.paths.qr_semiparam= fullfile(outDir, 'semi_param', 'semi_param_distr_inflation.mat');
cfg.paths.bvar        = fullfile(outDir, 'BVAR', 'BVAR_inflation_pred_q.mat');

%% ── Global settings ──────────────────────────────────────────────────────
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth', 1);
rng(0, 'twister');

addpath(fullfile(scriptDir, 'intermediate_codes'));
addpath(fullfile(scriptDir, 'functions'));
addpath(fullfile(scriptDir, 'functions', 'azzalini'));
addpath(fullfile(scriptDir, 'functions', 'CRPS'));

nHor    = numel(cfg.eval_horizons);
models  = {'fanchart', 'qr', 'qr_sp', 'bvar'};
nModels = numel(models);

% Shared deterministic quantile grid used by all three models for CRPS/wCRPS.
% G=999 points: accurate to O(1/G²) with no randomness.
tau_grid = (1:999)' / 1000;   % 999 × 1, values in (0,1)
% Right-tail weight kernel v(tau) = tau^2  [Gneiting-Ranjan QW=4]
v_right  = tau_grid.^2;

% Helper: uniform and right-tail wCRPS from a quantile grid q_grid and realisation y
%   CRPS(F,y)  = 2 * mean( (1{y<q} - tau) * (q - y) )
%   wCRPS(F,y) = 2 * mean( v(tau) * (1{y<q} - tau) * (q - y) )
crps_from_qgrid  = @(q, y) 2 * mean((double(y < q) - tau_grid) .* (q - y));
wcrps_from_qgrid = @(q, y) 2 * mean(v_right .* (double(y < q) - tau_grid) .* (q - y));

%% ════════════════════════════════════════════════════════════════════════
%%  SHARED ACTUAL DATA  (g4cpi — loaded once, aligned to common grid)
%% ════════════════════════════════════════════════════════════════════════
fprintf('Building shared actual data (g4cpi)...\n');

T_cpi      = readtable(cfg.paths.monthly, 'Sheet', 'g4cpi');
rawDates   = T_cpi{:, 1};
yearVec_c  = str2double(extractBetween(rawDates, 1, 4));
monthVec_c = str2double(extractAfter(rawDates, 'm'));
cpiDates   = datetime(yearVec_c, monthVec_c, 1);   % 1st of each month
ukColIdx   = find(strcmpi(T_cpi.Properties.VariableNames, 'UK'), 1);
cpiSeries  = T_cpi{:, ukColIdx};
cpi_ym     = year(cpiDates)*12 + month(cpiDates);  % integer year-month

% Common quarterly evaluation grid (quarter-end dates)
eval_origins = (cfg.eval_start : calmonths(3) : cfg.eval_end)';
nQOrig_eval  = numel(eval_origins);
eval_ym      = year(eval_origins)*12 + month(eval_origins);

% actualvar(ih, j) = g4cpi realisation cfg.eval_horizons(ih) quarters
% ahead of eval_origins(j).
actualvar = NaN(nHor, nQOrig_eval);
for j = 1:nQOrig_eval
    for ih = 1:nHor
        k      = cfg.eval_horizons(ih);
        tgt_ym = eval_ym(j) + 3*k;
        idx    = find(cpi_ym == tgt_ym, 1);
        if ~isempty(idx)
            actualvar(ih, j) = cpiSeries(idx);
        end
    end
end

clearvars T_cpi rawDates yearVec_c monthVec_c cpiDates ukColIdx cpiSeries;
% Keep: actualvar, eval_origins, nQOrig_eval, eval_ym, cpi_ym

%% ════════════════════════════════════════════════════════════════════════
%%  LOAD ALL MODELS AND COMPUTE PITs, CRPS, wCRPS, SHARPNESS
%%
%%  zinf_all{m}  — nHor × nQOrig_eval  PIT values
%%  crps_all{m}  — nHor × nQOrig_eval  uniform CRPS
%%  wcrps_all{m} — nHor × nQOrig_eval  right-tail weighted CRPS (v=tau^2)
%%  qntl_all{m}  — nHor × nQOrig_eval × nQL  quantiles (sharpness)
%%  originDT_all{m} = eval_origins  (identical across all models)
%% ════════════════════════════════════════════════════════════════════════

% Quantile levels used for sharpness interval widths
q_sharp_lev = [0.05, 0.10, 0.25, 0.35, 0.50, 0.65, 0.75, 0.90, 0.95];
nQL = numel(q_sharp_lev);   % 9 levels

zinf_all     = cell(nModels, 1);
crps_all     = cell(nModels, 1);
wcrps_all    = cell(nModels, 1);
qntl_all     = cell(nModels, 1);
originDT_all = cell(nModels, 1);

%% ── (1) Fan chart ────────────────────────────────────────────────────────
fprintf('Loading BOE fan chart...\n');

start_date = cfg.fanchart.start_date;
end_date   = cfg.fanchart.end_date;   %#ok<NASGU>
end_fcst   = cfg.fanchart.end_fcst;   %#ok<NASGU>

boefsctdata;   % produces mtestdata (13 × 3*nQOrig_fc)

% Extract parameter arrays (13 × nQOrig_fc)
modevar   = mtestdata(:, 1:3:end);
meanvar   = mtestdata(:, 2:3:end);
dispvar   = mtestdata(:, 3:3:end);
nQOrig_fc = size(meanvar, 2);

% Fan chart origin datetimes
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
        row = cfg.eval_horizons(ih) + 1;
        if row > size(meanvar, 1), continue; end
        if isnan(meanvar(row, fc_col)) || isnan(modevar(row, fc_col)), continue; end

        mn = meanvar(row, fc_col);
        mo = modevar(row, fc_col);
        dv = dispvar(row, fc_col);

        if mn ~= mo
            % Two-piece normal (asymmetric fan chart)
            gam = stdtogam(mn, mo, dv);
            sig = mom2g(dv, gam);
            [m_par, s1, s2] = momtopar(mn, mo, sig);

            % Quantile function of TPN:
            %   tau <= s1/(s1+s2): Q(tau) = m + s1*Phi^{-1}(tau*(s1+s2)/(2*s1))
            %   tau >  s1/(s1+s2): Q(tau) = m - s2*Phi^{-1}((1-tau)*(s1+s2)/(2*s2))
            p_left    = s1 / (s1 + s2);
            lm        = tau_grid <= p_left;
            q_grid    = NaN(size(tau_grid));
            q_grid( lm) = m_par + s1 .* norminv(tau_grid( lm) .* (s1+s2) ./ (2*s1));
            q_grid(~lm) = m_par - s2 .* norminv((1-tau_grid(~lm)) .* (s1+s2) ./ (2*s2));

            % Same formula at the specific sharpness quantile levels
            lm_s = q_sharp_lev <= p_left;
            q_s  = NaN(1, nQL);
            q_s( lm_s) = m_par + s1 .* norminv(q_sharp_lev( lm_s) .* (s1+s2) ./ (2*s1));
            q_s(~lm_s) = m_par - s2 .* norminv((1-q_sharp_lev(~lm_s)) .* (s1+s2) ./ (2*s2));

            % PIT via numerical integration of TPN PDF
            act = actualvar(ih, j);
            if ~isnan(act)
                zinf_fc(ih, j) = integral(@(y) ftp(y, m_par, s1, s2), -1e4, act, ...
                                          'RelTol', 1e-6, 'AbsTol', 1e-10);
                crps_fc(ih,j)  = crps_from_qgrid(q_grid, act);
                wcrps_fc(ih,j) = wcrps_from_qgrid(q_grid, act);
            end
        else
            % Symmetric (normal) fan chart
            q_grid = mn + dv .* norminv(tau_grid);
            q_s    = mn + dv .* norminv(q_sharp_lev);

            act = actualvar(ih, j);
            if ~isnan(act)
                zinf_fc(ih, j) = normcdf((act - mn) / dv);
                crps_fc(ih,j)  = crps_from_qgrid(q_grid, act);
                wcrps_fc(ih,j) = wcrps_from_qgrid(q_grid, act);
            end
        end

        % Sharpness: store quantiles (even when actual is NaN — widths don't need it)
        qntl_fc(ih, j, :) = q_s;
    end
end

zinf_all{1}     = zinf_fc;
crps_all{1}     = crps_fc;
wcrps_all{1}    = wcrps_fc;
qntl_all{1}     = qntl_fc;
originDT_all{1} = eval_origins;

clearvars modevar meanvar dispvar nQOrig_fc fc_origins fc_ym fc_col ...
          zinf_fc crps_fc wcrps_fc qntl_fc start_date end_date end_fcst ...
          gam sig ih j row m_par s1 s2 mn mo dv act p_left lm lm_s q_s q_grid mtestdata;

%% ── (2) RASS Skew-t ──────────────────────────────────────────────────────
fprintf('Loading RASS skew-t model...\n');

sd = load(cfg.paths.qr_sktparam, 'lc_skt', 'sc_skt', 'sh_skt', 'df_skt');
ad = load(cfg.paths.qr_actual,   'actual_var', 'idx_est', 'dateNumeric_full');

lc_qr = sd.lc_skt;
sc_qr = sd.sc_skt;
sh_qr = sd.sh_skt;
df_qr = sd.df_skt;

idx_est_qr = ad.idx_est;
dates_qr   = ad.dateNumeric_full;

nOrig_qr      = min(size(lc_qr, 1), numel(dates_qr) - idx_est_qr + 1);
orig_dates_qr = datetime(dates_qr(idx_est_qr : idx_est_qr + nOrig_qr - 1), ...
                          'ConvertFrom', 'datenum');
qr_ym         = year(orig_dates_qr)*12 + month(orig_dates_qr);

zinf_qr  = NaN(nHor, nQOrig_eval);
crps_qr  = NaN(nHor, nQOrig_eval);
wcrps_qr = NaN(nHor, nQOrig_eval);
qntl_qr  = NaN(nHor, nQOrig_eval, nQL);

for j = 1:nQOrig_eval
    t_idx = find(qr_ym == eval_ym(j), 1);
    if isempty(t_idx), continue; end

    for ih = 1:nHor
        k         = cfg.eval_horizons(ih);
        h_monthly = 3*k + 1;
        if h_monthly > size(lc_qr, 2), continue; end

        lc_v = lc_qr(t_idx, h_monthly);
        sc_v = sc_qr(t_idx, h_monthly);
        sh_v = sh_qr(t_idx, h_monthly);
        df_v = df_qr(t_idx, h_monthly);
        if any(isnan([lc_v, sc_v, sh_v, df_v])) || sc_v <= 0, continue; end

        % Build CDF on an x-grid wide enough for the right tail of the skew-t.
        % This avoids qskt/pskt which both fail for large shape parameters.
        nG   = 2000;
        x_lo = lc_v - 15 * sc_v;
        x_hi = lc_v + max(15, 8 * abs(sh_v)) * sc_v;
        xg   = linspace(x_lo, x_hi, nG)';
        pdg  = dskt(xg, lc_v, sc_v, sh_v, df_v);
        cdg  = min(max(cumtrapz(xg, pdg), 0), 1);   % CDF, clamped to [0,1]

        % Sharpness: interpolate quantile levels from CDF grid (no actual needed)
        [cu, ia] = unique(cdg, 'stable');
        xu = xg(ia);
        qntl_qr(ih, j, :) = interp1(cu, xu, q_sharp_lev, 'linear', NaN);

        % PIT, CRPS, wCRPS only when actual is available
        act = actualvar(ih, j);
        if isnan(act), continue; end

        % PIT via numerical integration of skew-t PDF
        zinf_qr(ih, j) = integral( ...
            @(x) dskt(x, lc_v, sc_v, sh_v, df_v), -1e4, act, ...
            'RelTol', 1e-6, 'AbsTol', 1e-10);

        % CRPS  = integral (F(x) - 1{x>=y})^2 dx
        crps_qr(ih,j)  = trapz(xg, (cdg - double(xg >= act)).^2);
        % wCRPS right tail: v(tau)=tau^2, substituting tau=F(x):
        %   wCRPS = 2 * integral F(x)^2 * (1{x>y} - F(x)) * (x-y) * f(x) dx
        wcrps_qr(ih,j) = 2 * trapz(xg, ...
            cdg.^2 .* (double(xg > act) - cdg) .* (xg - act) .* pdg);
    end
end

zinf_all{2}     = zinf_qr;
crps_all{2}     = crps_qr;
wcrps_all{2}    = wcrps_qr;
qntl_all{2}     = qntl_qr;
originDT_all{2} = eval_origins;

clearvars sd ad lc_qr sc_qr sh_qr df_qr idx_est_qr dates_qr ...
          nOrig_qr orig_dates_qr qr_ym zinf_qr crps_qr wcrps_qr qntl_qr ...
          ih j k h_monthly lc_v sc_v sh_v df_v t_idx act ...
          nG x_lo x_hi xg pdg cdg;

%% ── (3) RASS Semi-parametric (Mitchell-Poon-Zhu) ────────────────────────
%  semi_param_distr: nOrigins × 20000 × 37  (Monte Carlo draws; only
%  quarterly-origin rows are populated, rest are NaN).
fprintf('Loading RASS semi-parametric model...\n');

sp = load(cfg.paths.qr_semiparam, 'semi_param_distr');
semi_param_distr = sp.semi_param_distr;   % nOrigins × 20000 × 37

% Reconstruct monthly origin dates from auxiliary file.
% dateNumeric_full covers startT→endT; origins start at idx_est.
% Take from idx_est to end, then clip to however many rows the matrix has.
ad_sp     = load(cfg.paths.qr_actual, 'idx_est', 'dateNumeric_full');
full_dn   = ad_sp.dateNumeric_full(ad_sp.idx_est : end);
n_sp_use  = min(numel(full_dn), size(semi_param_distr, 1));
sp_orig_dates = datetime(full_dn(1:n_sp_use), 'ConvertFrom', 'datenum');
sp_ym = year(sp_orig_dates)*12 + month(sp_orig_dates);

zinf_sp  = NaN(nHor, nQOrig_eval);
crps_sp  = NaN(nHor, nQOrig_eval);
wcrps_sp = NaN(nHor, nQOrig_eval);
qntl_sp  = NaN(nHor, nQOrig_eval, nQL);

for j = 1:nQOrig_eval
    t_idx = find(sp_ym == eval_ym(j), 1);
    if isempty(t_idx), continue; end

    for ih = 1:nHor
        k         = cfg.eval_horizons(ih);
        h_monthly = 3*k + 1;
        if h_monthly > size(semi_param_distr, 3), continue; end

        draws = squeeze(semi_param_distr(t_idx, :, h_monthly))';   % 20000 × 1
        draws = draws(~isnan(draws));
        if numel(draws) < 100, continue; end

        % Sharpness: quantiles from empirical draw distribution
        qntl_sp(ih, j, :) = quantile(draws, q_sharp_lev);

        % PIT, CRPS, wCRPS only when actual is available
        act = actualvar(ih, j);
        if isnan(act), continue; end

        % PIT = empirical CDF at actual value (fraction of draws ≤ actual)
        zinf_sp(ih, j) = mean(draws <= act);

        % Quantile grid at 999 points from empirical distribution
        q_grid = quantile(draws, tau_grid);

        crps_sp(ih, j)  = crps_from_qgrid(q_grid, act);
        wcrps_sp(ih, j) = wcrps_from_qgrid(q_grid, act);
    end
end

zinf_all{3}     = zinf_sp;
crps_all{3}     = crps_sp;
wcrps_all{3}    = wcrps_sp;
qntl_all{3}     = qntl_sp;
originDT_all{3} = eval_origins;

clearvars sp semi_param_distr ad_sp full_dn n_sp_use sp_orig_dates sp_ym ...
          zinf_sp crps_sp wcrps_sp qntl_sp ...
          ih j k h_monthly draws act t_idx q_grid;

%% ── (4) BVAR  (Normal predictive — mean and std from Gibbs draws) ────────
%  The BVAR posterior predictive at each horizon is characterised by
%  pred_mu (posterior mean of draws) and pred_sigma (posterior std dev).
%  These are stored directly in BVAR_inflation.m, so no quantile-fitting
%  is needed here.  PIT and CRPS use the Normal distribution analytically.
fprintf('Loading BVAR...\n');

bd = load(cfg.paths.bvar, 'pred_mu', 'pred_sigma', 'dateNumeric_est');
pred_mu_bv    = bd.pred_mu;      % nOrigins × 37
pred_sigma_bv = bd.pred_sigma;   % nOrigins × 37

orig_dates_bv = datetime(bd.dateNumeric_est, 'ConvertFrom', 'datenum')';
bv_ym         = year(orig_dates_bv)*12 + month(orig_dates_bv);

zinf_bv  = NaN(nHor, nQOrig_eval);
crps_bv  = NaN(nHor, nQOrig_eval);
wcrps_bv = NaN(nHor, nQOrig_eval);
qntl_bv  = NaN(nHor, nQOrig_eval, nQL);

for j = 1:nQOrig_eval
    t_idx = find(bv_ym == eval_ym(j), 1);
    if isempty(t_idx), continue; end

    for ih = 1:nHor
        k         = cfg.eval_horizons(ih);
        h_monthly = 3*k + 1;
        if h_monthly > size(pred_mu_bv, 2), continue; end

        mu_bv    = pred_mu_bv(t_idx, h_monthly);
        sigma_bv = pred_sigma_bv(t_idx, h_monthly);
        if isnan(mu_bv) || isnan(sigma_bv) || sigma_bv <= 0, continue; end

        % Sharpness: Normal quantiles at specific levels (no actual needed)
        qntl_bv(ih, j, :) = norminv(q_sharp_lev, mu_bv, sigma_bv);

        % PIT, CRPS, wCRPS only when actual is available
        act = actualvar(ih, j);
        if isnan(act), continue; end

        % PIT via Normal CDF
        zinf_bv(ih, j) = normcdf(act, mu_bv, sigma_bv);

        % Quantile grid via Normal inverse CDF
        q_grid = norminv(tau_grid, mu_bv, sigma_bv);

        crps_bv(ih,j)  = crps_from_qgrid(q_grid, act);
        wcrps_bv(ih,j) = wcrps_from_qgrid(q_grid, act);
    end
end

zinf_all{4}     = zinf_bv;
crps_all{4}     = crps_bv;
wcrps_all{4}    = wcrps_bv;
qntl_all{4}     = qntl_bv;
originDT_all{4} = eval_origins;

clearvars bd pred_mu_bv pred_sigma_bv orig_dates_bv bv_ym ...
          zinf_bv crps_bv wcrps_bv qntl_bv ih j k h_monthly act ...
          mu_bv sigma_bv t_idx q_grid;

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
    results.(models{m}).wcrps_mean   = NaN(nHor, 1);   % right-tail weighted CRPS  [v(tau)=tau^2]
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

        % wCRPS mean (right-tail weighted, v=tau^2)
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
% MC     = 1000;
% bootMC = 1000;
% rng(bootMC, 'twister');
% rvec_gmm = 0:0.001:1;
% 
% for m = 1:nModels
%     mod  = models{m};
%     zinf = zinf_all{m};
%     z    = zinf';
%     P    = size(z, 1);
%     el   = floor(P^(1/4));
%     Hz   = size(z, 2);
%     KS   = results.(mod).KS_vec(:)';
%     CVM  = results.(mod).CVM_vec(:)';
% 
%     QVrejvecs       = zeros(MC, 3);
%     CVMrejvecs      = zeros(MC, 3);
%     QVrejvecs_bonf  = zeros(MC, 3);
%     CVMrejvecs_bonf = zeros(MC, 3);
% 
%     parfor j = 1:MC
%         stream1 = RandStream('mrg32k3a', 'seed', 4829575);
%         stream1.Substream = j;
%         [QVrej_j, CVMrej_j, QVbonf_j, CVMbonf_j] = ...
%             size_statistic_h2(z, KS, CVM, Hz, stream1, rvec_gmm, el, bootMC);
%         QVrejvecs(j,:)       = QVrej_j;
%         CVMrejvecs(j,:)      = CVMrej_j;
%         QVrejvecs_bonf(j,:)  = QVbonf_j;
%         CVMrejvecs_bonf(j,:) = CVMbonf_j;
%     end
% 
%     results.(mod).gmm_ks       = mean(QVrejvecs,       1);
%     results.(mod).gmm_cvm      = mean(CVMrejvecs,       1);
%     results.(mod).gmm_ks_bonf  = mean(QVrejvecs_bonf,  1);
%     results.(mod).gmm_cvm_bonf = mean(CVMrejvecs_bonf, 1);
% 
%     fprintf('GMM done: %s\n', mod);
% end

%% ════════════════════════════════════════════════════════════════════════
%%  SAVE COMPARISON RESULTS
%% ════════════════════════════════════════════════════════════════════════

horiz    = cfg.eval_horizons(:);
filename = fullfile(evalDir, 'comparison_fcst_eval_inflation.xlsx');

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
writetable(buildTab('wcrps_mean'),    filename, 'Sheet', 'wCRPS_right_mean');

% gmm_ks_mat  = cell2mat(cellfun(@(m) results.(m).gmm_ks(:)',  models(:)', 'UniformOutput', false)')';
% gmm_cvm_mat = cell2mat(cellfun(@(m) results.(m).gmm_cvm(:)', models(:)', 'UniformOutput', false)')';
% gmm_tab = table(models', gmm_ks_mat, gmm_cvm_mat, ...
%     'VariableNames', {'model','GMM_KS_std_bonf_w_wi','GMM_CVM_std_bonf_w_wi'});
% writetable(gmm_tab, filename, 'Sheet', 'GMM_summary');

%% ════════════════════════════════════════════════════════════════════════
%%  PRINT SUMMARY TABLES  (one per model)
%% ════════════════════════════════════════════════════════════════════════

model_labels = {'Fan Chart (BOE)', 'QR Skew-t', 'QR Semi-param', 'BVAR'};

sel_h   = [0 1 2 3 4 8 12];
sel_idx = arrayfun(@(h) find(cfg.eval_horizons == h, 1), sel_h);

frcode = @(v) char(70*(~isnan(v) & v < 0.5) + ...
                   82*(~isnan(v) & v >= 0.5) + ...
                   45*isnan(v));

hdr = ['  h   | Knu pval  |  MW pval  | CRPS mean | wCRPS(R)' ...
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

pstar = @(p) repmat('*', 1, (p<0.10) + (p<0.05) + (p<0.01));

col_spec = ['l', repmat('c', 1, nModels)];
col_hdrs = strjoin(cellfun(@(l) ['\textbf{', strrep(l,' ','\;'), '}'], ...
                   model_labels, 'UniformOutput', false), ' & ');

latex_file_infl = fullfile(evalDir, 'table_fcst_eval_inflation.tex');
fid = fopen(latex_file_infl, 'w');

fprintf(fid, '%% Inflation forecast evaluation table — auto-generated\n');
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Inflation Forecast Evaluation}\n');
fprintf(fid, '\\label{tab:fcst_eval_inflation}\n');
fprintf(fid, '\\small\n');
fprintf(fid, '\\begin{tabular}{%s}\n', col_spec);
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Horizon & %s \\\\\n', col_hdrs);
fprintf(fid, '\\midrule\n');

tex_blocks_infl = { ...
    'Knueppel (2015) p-value', 'knueppel_pval', true;  ...
    'Mitchell-Weale p-value',  'MW_pval',       true;  ...
    'CRPS (mean)',             'crps_mean',     false; ...
    'wCRPS right tail (mean)', 'wcrps_mean',    false  ...
};

for tb = 1:size(tex_blocks_infl, 1)
    blk_label  = tex_blocks_infl{tb, 1};
    blk_field  = tex_blocks_infl{tb, 2};
    blk_ispval = tex_blocks_infl{tb, 3};
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
fprintf(fid, '\\item Notes: $h$ is horizon in months. ');
fprintf(fid, 'p-values: $^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$. ');
fprintf(fid, 'wCRPS uses right-tail weight $v(\\tau)=\\tau^2$. ');
fprintf(fid, 'Covid origins excluded.\n');
fprintf(fid, '\\end{tablenotes}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);

% Print to console
fprintf('\n\n%% ── LaTeX table (inflation) ─────────────────────────────\n');
fid2 = fopen(latex_file_infl, 'r');
while ~feof(fid2), fprintf('%s\n', fgetl(fid2)); end
fclose(fid2);
fprintf('\nLaTeX table saved to %s\n', latex_file_infl);

%% ════════════════════════════════════════════════════════════════════════
%%  SHARPNESS DIAGNOSTIC
%%
%%  For each model and horizon, compute interval widths across forecast
%%  origins, then plot median ± 1 SD by horizon (following the style of
%%  sharpness_test_inflation_v2.m but for all three models together).
%%
%%  Quantile indices in qntl_all{m}(ih, j, :):
%%    1=q05  2=q10  3=q25  4=q35  5=q50  6=q65  7=q75  8=q90  9=q95
%% ════════════════════════════════════════════════════════════════════════
fprintf('\nComputing sharpness widths...\n');

h_vec_sharp = cfg.eval_horizons(:);   % 0..12 quarters

% Width definitions (indices into q_sharp_lev = [.05 .10 .25 .35 .50 .65 .75 .90 .95])
%   w50  = q75 - q25  → cols 7 - 3
%   w90  = q95 - q05  → cols 9 - 1
%   w30  = q65 - q35  → cols 6 - 4
%   wb20 = q25 - q05  → cols 3 - 1  (bottom 20%)
%   wt20 = q95 - q75  → cols 9 - 7  (top 20%)

% 2×2 layout: Central 30% | Central 50% | Bottom 20% | Top 20%
sharp_labels = {'Central 30%','Central 50%','Bottom 20%','Top 20%'};
nW = numel(sharp_labels);

% mu_w{m}(ih, iw) = median width across origins at horizon ih, width type iw
mu_w  = cell(nModels, 1);
sd_w  = cell(nModels, 1);

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

% Print sharpness table (median widths at selected horizons)
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

% ── Shared style (Fan Chart, QR Skew-t, BVAR only in charts) ─────────────
chart_m   = {'fanchart','qr','bvar'};   % subset for plotting
chart_lbl = {'Fan Chart (BOE)','QR Skew-t','BVAR'};
nCM       = numel(chart_m);
colors_sh  = [0.12 0.47 0.71; 0.84 0.15 0.16; 0.17 0.63 0.17];
fills_sh   = [0.70 0.80 1.00; 1.00 0.80 0.80; 0.80 1.00 0.80];
lssh       = {'-', '--', ':'};
markers_sh = {'o', 's', '^'};
lwsh       = 1.8;
h_vec      = cfg.eval_horizons(:);

% Map chart_m → index in mu_w / model_labels
chart_idx = cellfun(@(c) find(strcmp(models, c)), chart_m);

%% ── CHART 1: Sharpness — 2×2 ─────────────────────────────────────────────
fig_sh = figure('Name','Sharpness — CPI Inflation','NumberTitle','off', ...
                'Color','w','Position',[120 120 900 750]);
tl_sh  = tiledlayout(fig_sh, 2, 2, 'TileSpacing','compact','Padding','loose');
sgtitle(fig_sh,'CPI Inflation — Forecast Sharpness (median width ± 1 SD)', ...
        'FontSize',12,'FontWeight','bold');

legH_sh = gobjects(nCM,1);
for iw = 1:nW
    ax = nexttile(tl_sh, iw);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for ci = 1:nCM
        m    = chart_idx(ci);
        mu_v = mu_w{m}(:, iw);
        sd_v = sd_w{m}(:, iw);
        xx   = [h_vec_sharp; flipud(h_vec_sharp)];
        yy   = [mu_v - sd_v; flipud(mu_v + sd_v)];
        patch('Parent',ax,'XData',xx,'YData',yy, ...
              'FaceColor',fills_sh(ci,:),'EdgeColor','none', ...
              'FaceAlpha',0.40,'HandleVisibility','off');
        hh = plot(ax, h_vec_sharp, mu_v, lssh{ci}, ...
                  'Color',colors_sh(ci,:),'LineWidth',lwsh, ...
                  'Marker',markers_sh{ci},'MarkerSize',4, ...
                  'MarkerFaceColor',colors_sh(ci,:), ...
                  'DisplayName',chart_lbl{ci});
        if iw == 1, legH_sh(ci) = hh; end
    end
    xlabel(ax,'Horizon (months ahead)','FontSize',9);
    ylabel(ax,'Width (pp)','FontSize',9);
    title(ax, sharp_labels{iw},'FontSize',10,'FontWeight','bold');
    set(ax,'XTick',h_vec_sharp,'Box','on','FontSize',9);
    hold(ax,'off');
end
lgd_sh = legend(legH_sh, chart_lbl,'Orientation','horizontal','NumColumns',3);
try, lgd_sh.Layout.Tile = 'south'; catch
    set(lgd_sh,'Position',[0.25 0.01 0.50 0.04]); end
exportgraphics(fig_sh, fullfile(evalDir,'sharpness_fcst_eval_inflation.png'),'Resolution',300);
fprintf('\nSharpness chart saved.\n');

%% ── CHART 2: Calibration p-values — 1×3 ─────────────────────────────────
fig_pv = figure('Name','Calibration Tests — CPI Inflation','NumberTitle','off', ...
                'Color','w','Position',[80 80 900 420]);
tl_pv  = tiledlayout(fig_pv, 1, 2,'TileSpacing','compact','Padding','loose');
sgtitle(fig_pv,'CPI Inflation — Calibration Tests (p-values)', ...
        'FontSize',12,'FontWeight','bold');

pval_fields = {'knueppel_pval','MW_pval'};
pval_titles = {'Knueppel (2015) p-value','Mitchell-Weale (2023) p-value'};
legH_pv = gobjects(nCM,1);

for t = 1:2
    ax = nexttile(tl_pv, t);
    hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for ci = 1:nCM
        hh = plot(ax, h_vec, results.(chart_m{ci}).(pval_fields{t}), ...
             'Color',colors_sh(ci,:),'LineStyle',lssh{ci},'LineWidth',lwsh, ...
             'Marker',markers_sh{ci},'MarkerSize',4,'MarkerFaceColor',colors_sh(ci,:), ...
             'DisplayName',chart_lbl{ci});
        if t == 1, legH_pv(ci) = hh; end
    end
    yline(ax,0.10,'Color',[0.4 0.4 0.4],'LineStyle',':','LineWidth',1.0, ...
          'Label','10%','LabelVerticalAlignment','bottom','HandleVisibility','off');
    yline(ax,0.05,'Color',[0.2 0.2 0.2],'LineStyle','--','LineWidth',1.0, ...
          'Label','5%','LabelVerticalAlignment','bottom','HandleVisibility','off');
    ylim(ax,[0 1]);
    xlabel(ax,'Horizon (months ahead)','FontSize',9);
    title(ax, pval_titles{t},'FontSize',10,'FontWeight','bold');
    set(ax,'XTick',h_vec,'Box','on','FontSize',9);
    hold(ax,'off');
end
lgd_pv = legend(legH_pv, chart_lbl,'Orientation','horizontal','NumColumns',3);
try, lgd_pv.Layout.Tile = 'south'; catch
    set(lgd_pv,'Position',[0.20 0.01 0.60 0.04]); end
exportgraphics(fig_pv, fullfile(evalDir,'pvalues_fcst_eval_inflation.png'),'Resolution',300);
fprintf('\nCalibration p-value chart saved.\n');

%% ── CHART 3: Cumulative CRPS & wCRPS — 2×2 ──────────────────────────────
cum_h       = [1, 4];
cum_h_names = {'1M ahead (h=1)', '1Y ahead (h=12)'};

fig_sc = figure('Name','Cumulative Scores — CPI Inflation','NumberTitle','off', ...
                'Color','w','Position',[100 100 1000 750]);
tl_sc  = tiledlayout(fig_sc, 2, 2,'TileSpacing','compact','Padding','loose');
sgtitle(fig_sc,'CPI Inflation — Cumulative CRPS & wCRPS (right tail)', ...
        'FontSize',12,'FontWeight','bold');

legH_sc = gobjects(nCM,1);

% Row 1: Cumulative CRPS at h=1 and h=12
for ip = 1:2
    ax_c = nexttile(tl_sc, ip);
    hold(ax_c,'on'); box(ax_c,'on'); grid(ax_c,'on');
    ih_row = find(cfg.eval_horizons == cum_h(ip), 1);
    for ci = 1:nCM
        score_row = crps_all{chart_idx(ci)}(ih_row, :);
        valid = ~isnan(score_row);
        hh = plot(ax_c, eval_origins, cumsum(score_row .* valid), ...
             'Color',colors_sh(ci,:),'LineStyle',lssh{ci},'LineWidth',lwsh, ...
             'DisplayName',chart_lbl{ci});
        if ip == 1, legH_sc(ci) = hh; end
    end
    xlabel(ax_c,'Evaluation origin','FontSize',9);
    ylabel(ax_c,'Cumulative CRPS','FontSize',9);
    title(ax_c, sprintf('Cum. CRPS — %s', cum_h_names{ip}), ...
          'FontSize',10,'FontWeight','bold');
    set(ax_c,'FontSize',9);
    datetick(ax_c,'x','yyyy','keepticks','keeplimits'); xtickangle(ax_c,45);
    hold(ax_c,'off');
end

% Row 2: Cumulative wCRPS (right tail) at h=1 and h=12
for ip = 1:2
    ax_w = nexttile(tl_sc, ip+2);
    hold(ax_w,'on'); box(ax_w,'on'); grid(ax_w,'on');
    ih_row = find(cfg.eval_horizons == cum_h(ip), 1);
    for ci = 1:nCM
        score_row = wcrps_all{chart_idx(ci)}(ih_row, :);
        valid = ~isnan(score_row);
        plot(ax_w, eval_origins, cumsum(score_row .* valid), ...
             'Color',colors_sh(ci,:),'LineStyle',lssh{ci},'LineWidth',lwsh, ...
             'DisplayName',chart_lbl{ci});
    end
    xlabel(ax_w,'Evaluation origin','FontSize',9);
    ylabel(ax_w,'Cumulative wCRPS','FontSize',9);
    title(ax_w, sprintf('Cum. wCRPS (right tail) — %s', cum_h_names{ip}), ...
          'FontSize',10,'FontWeight','bold');
    set(ax_w,'FontSize',9);
    datetick(ax_w,'x','yyyy','keepticks','keeplimits'); xtickangle(ax_w,45);
    hold(ax_w,'off');
end

lgd_sc = legend(legH_sc, chart_lbl,'Orientation','horizontal','NumColumns',3);
try, lgd_sc.Layout.Tile = 'south'; catch
    set(lgd_sc,'Position',[0.20 0.01 0.60 0.04]); end
exportgraphics(fig_sc, fullfile(evalDir,'scores_fcst_eval_inflation.png'),'Resolution',300);
fprintf('\nCumulative scores chart saved.\n');

fprintf('\n── DONE: comparison saved to %s ──\n', filename);
