%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FORECAST EVALUATION — CPI INFLATION  (multi-model comparison)
%%
%%  Authors: David Aikman, Rhys Bidder, Simon Lloyd, Giulia Mantoan,
%%           Simone Maso, Aditya Mori, Matthew Tong
%%
%%  Tests: KS · Rossi-Sekhposyan (2019) · Berkowitz (2001) ·
%%         Knueppel (2015) · Mitchell-Weale (2023) · Galvao-Mantoan-Mitchell
%%
%%  Compares: Fan Chart (BOE) · QR (RASS) · BVAR
%%  Results are saved side-by-side for each test statistic.
%%
%%  Notes:
%%    - Evaluation is on QUARTERLY origins (months 3,6,9,12) and
%%      QUARTERLY horizons (0,3,6,...,36 months ahead = 0..12 quarters ahead).
%%    - No Covid exclusion applied to inflation.
%%    - Fan chart actual data: g4cpi from IaRDataRaw_monthly_M.xlsx,
%%      filtered to quarter-end months; rows 1:3:37 select quarterly horizons.
%%    - QR/BVAR monthly origins filtered to quarter-ends; monthly horizon
%%      h_monthly = 3*k+1 maps to quarterly horizon k (0-indexed).
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

% Common evaluation window (set [] to use each model's full available sample)
cfg.eval_start = datetime('31-Mar-2010', 'InputFormat', 'dd-MMM-yyyy');
cfg.eval_end   = datetime('30-Sep-2022', 'InputFormat', 'dd-MMM-yyyy');

% No Covid exclusion for inflation
cfg.covid_exclude = false;

% Fan chart settings
cfg.fanchart.start_date = datetime('31-Mar-2010', 'InputFormat', 'dd-MMM-yyyy');
cfg.fanchart.end_date   = datetime('30-Sep-2022', 'InputFormat', 'dd-MMM-yyyy');
cfg.fanchart.end_fcst   = datetime('30-Sep-2025', 'InputFormat', 'dd-MMM-yyyy');

% Root directory that boefsctdata uses via cd() to locate
% Data\cpi_infl_projection_parameters_mpc.xlsx
cfg.paths.boe_root  = scriptDir;
cfg.paths.monthly   = fullfile(scriptDir, 'IaRDataRaw_monthly_M.xlsx');
cfg.paths.qr_pred_q = fullfile(outDir, 'predicted_quantiles_inflation_OOS.mat');
cfg.paths.qr_actual = fullfile(outDir, 'actual_inflation_mom_OOS.mat');
cfg.paths.bvar      = fullfile(outDir, 'BVAR', 'BVAR_inflation_pred_q.mat');

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

%% ════════════════════════════════════════════════════════════════════════
%%  LOAD ALL MODELS AND COMPUTE PITs
%%
%%  zinf_all{m} — nHor × nOrigins_m  matrix of PITs for model m.
%%  Row ih = cfg.eval_horizons(ih) quarters ahead.
%%  Columns not in cfg.eval_start/end window are left as NaN.
%% ════════════════════════════════════════════════════════════════════════

zinf_all     = cell(nModels, 1);
crps_all     = cell(nModels, 1);
originDT_all = cell(nModels, 1);

%% ── (1) Fan chart ────────────────────────────────────────────────────────
fprintf('Loading BOE fan chart...\n');

% boefsctdata has no clearvars — safe to call without workspace backup.
% It requires: start_date, end_date, end_fcst (datetimes) and cd() pointing
% to the folder that contains the Data\ subfolder.
start_date = cfg.fanchart.start_date;
end_date   = cfg.fanchart.end_date;   %#ok<NASGU>  (used inside boefsctdata)
end_fcst   = cfg.fanchart.end_fcst;

cd(cfg.paths.boe_root);
boefsctdata;   % produces mtestdata (13 × 3*nQOrig_fc); no clearvars

% ── Load actual CPI data from monthly file ────────────────────────────────
% Parse g4cpi sheet (date format: '1970m1')
T_cpi    = readtable(cfg.paths.monthly, 'Sheet', 'g4cpi');
rawDates = T_cpi{:, 1};
yearVec  = str2double(extractBetween(rawDates, 1, 4));
monthVec = str2double(extractAfter(rawDates, 'm'));
cpiDates = datetime(yearVec, monthVec, 1);   % first of each month

ukColIdx = find(strcmpi(T_cpi.Properties.VariableNames, 'UK'), 1);
cpiData  = T_cpi{:, ukColIdx};

% Filter to [start_date, end_fcst] at month level
ym_data  = year(cpiDates)*12 + month(cpiDates);
ym_start = year(start_date)*12 + month(start_date);
ym_fcst  = year(end_fcst)*12  + month(end_fcst);
keepIdx  = (ym_data >= ym_start) & (ym_data <= ym_fcst);
cpiData  = cpiData(keepIdx);
cpiDates = cpiDates(keepIdx);

% Build lag matrix (37 × n_monthly): row h = (h-1) months ahead
n_monthly = numel(cpiData);
interm    = NaN(n_monthly, n_monthly);
for ii = 1:n_monthly
    interm(1:(n_monthly-ii+1), ii) = cpiData(ii:end);
end
actualvar_monthly = interm(1:37, :);   % 37 × n_monthly
clear interm ii;

% Filter to quarter-end origins (months 3,6,9,12)
% and quarterly horizons (rows 1,4,7,...,37 → 0,1,...,12 quarters ahead)
qtr_mask_fc = ismember(month(cpiDates), [3 6 9 12]);
actualvar   = actualvar_monthly(1:3:37, qtr_mask_fc);   % 13 × nQOrig_cpi
cpiDates_qe = dateshift(cpiDates(qtr_mask_fc), 'end', 'month');   % end-of-month

% Fan chart PITs (Two-piece normal / Normal)
modevar  = mtestdata(:, 1:3:end);   % 13 × nQOrig_fc
meanvar  = mtestdata(:, 2:3:end);
dispvar  = mtestdata(:, 3:3:end);
nOrig_fc = size(meanvar, 2);

% Origin datetimes: fan chart quarters start from start_date
originDT_all{1} = (start_date + calmonths(3*(0:nOrig_fc-1)))';

zinf_fc = NaN(nHor, nOrig_fc);
crps_fc = NaN(nHor, nOrig_fc);
nSamp   = 20000;
for ih = 1:nHor
    row = cfg.eval_horizons(ih) + 1;   % 1-indexed row in meanvar / actualvar
    if row > size(meanvar, 1), continue; end
    for j = 1:nOrig_fc
        if j > size(actualvar, 2), continue; end
        if isnan(meanvar(row,j)) || isnan(modevar(row,j)), continue; end
        act = actualvar(row,j);
        if meanvar(row,j) ~= modevar(row,j)
            gam = stdtogam(meanvar(row,j), modevar(row,j), dispvar(row,j));
            sig = mom2g(dispvar(row,j), gam);
            [m_par, s1, s2] = momtopar(meanvar(row,j), modevar(row,j), sig);
            if ~isnan(act)
                zinf_fc(ih,j) = integral(@(y) ftp(y,m_par,s1,s2), -3, act);
                % CRPS: sample from two-piece normal (split-normal)
                u_samp = rand(nSamp, 1) < s1/(s1+s2);
                samp   = NaN(nSamp, 1);
                samp( u_samp) = m_par - abs(randn(sum( u_samp), 1)) .* s1;
                samp(~u_samp) = m_par + abs(randn(sum(~u_samp), 1)) .* s2;
                crps_fc(ih,j) = crps(samp, act, 1);
            end
        else
            if ~isnan(act)
                zinf_fc(ih,j) = normcdf((act - meanvar(row,j)) / dispvar(row,j));
                % CRPS: sample from symmetric normal
                samp = meanvar(row,j) + dispvar(row,j) * randn(nSamp, 1);
                crps_fc(ih,j) = crps(samp, act, 1);
            end
        end
    end
end
keep_fc_cols = ~all(isnan(zinf_fc), 1);
zinf_fc = zinf_fc(:, keep_fc_cols);
crps_fc = crps_fc(:, keep_fc_cols);
zinf_all{1} = zinf_fc;
crps_all{1} = crps_fc;

clearvars modevar meanvar dispvar zinf_fc crps_fc nSamp actualvar actualvar_monthly mtestdata ...
          gam sig ih row j m_par s1 s2 act u_samp samp keep_fc_cols interm n_monthly ...
          keepIdx nOrig_fc T_cpi rawDates yearVec monthVec cpiData cpiDates cpiDates_qe ...
          qtr_mask_fc ym_data ym_start ym_fcst ukColIdx;

%% ── (2) QR ───────────────────────────────────────────────────────────────
fprintf('Loading QR model...\n');

qd = load(cfg.paths.qr_pred_q, 'pred_q');
ad = load(cfg.paths.qr_actual,  'actual_var', 'idx_est', 'dateNumeric_full');
pred_q_qr  = qd.pred_q;       % nOrigins × Qn × 37 monthly horizons
actual_qr  = ad.actual_var;   % full monthly time series vector
idx_est_qr = ad.idx_est;
dates_qr   = ad.dateNumeric_full;
quant_qr   = (0.05:0.05:0.95)';

% Filter QR monthly origins to quarter-ends
% Clamp nOrig_qr to data actually present in dates_qr (pred_q may be
% pre-allocated larger than the data via date-arithmetic nOrigins)
nOrig_qr      = min(size(pred_q_qr, 1), numel(dates_qr) - idx_est_qr + 1);
orig_dates_qr = datetime(dates_qr(idx_est_qr : idx_est_qr + nOrig_qr - 1), ...
                          'ConvertFrom', 'datenum');
qtr_mask_qr   = ismember(month(orig_dates_qr), [3 6 9 12]);
qtr_idx_qr    = find(qtr_mask_qr);
nQOrig_qr     = numel(qtr_idx_qr);
originDT_all{2} = dateshift(orig_dates_qr(qtr_mask_qr)', 'end', 'month');

% PITs: for quarterly horizon k, monthly pred_q index h_monthly = 3*k+1
zinf_qr = NaN(nHor, nQOrig_qr);
crps_qr = NaN(nHor, nQOrig_qr);
for qi = 1:nQOrig_qr
    t_idx = qtr_idx_qr(qi);
    for ih = 1:nHor
        k          = cfg.eval_horizons(ih);        % quarters ahead (0-indexed)
        h_monthly  = 3*k + 1;                      % monthly index in pred_q
        if h_monthly > size(pred_q_qr, 3), continue; end
        q_vals = squeeze(pred_q_qr(t_idx, :, h_monthly));
        if all(isnan(q_vals)), continue; end
        tgt = (idx_est_qr + t_idx - 1) + (h_monthly - 1);
        if tgt > length(actual_qr), continue; end
        actual_h = actual_qr(tgt);
        if isnan(actual_h), continue; end
        pit = interp1(q_vals(:), quant_qr, actual_h, 'linear', 'extrap');
        zinf_qr(ih, qi) = min(max(pit, 0), 1);
        % CRPS: generate samples from smoothed quantile distribution
        samp_qr = QR_sm(q_vals(:)', quant_qr(:)');
        crps_qr(ih, qi) = crps(samp_qr(:), actual_h, 1);
    end
end
zinf_all{2} = zinf_qr;
crps_all{2} = crps_qr;

clearvars qd ad pred_q_qr actual_qr idx_est_qr dates_qr quant_qr zinf_qr crps_qr samp_qr ...
          ih qi k h_monthly q_vals tgt actual_h pit qtr_mask_qr qtr_idx_qr ...
          orig_dates_qr nOrig_qr nQOrig_qr;

%% ── (3) BVAR ─────────────────────────────────────────────────────────────
fprintf('Loading BVAR...\n');

bd = load(cfg.paths.bvar, 'pred_q', 'actualvar', 'dateNumeric_est', 'cfg');
pred_q_bv = bd.pred_q;       % nOrigins × Qn × 37 monthly horizons
actual_bv = bd.actualvar;    % nOrigins × 37  (already aligned)
quant_bv  = bd.cfg.quantiles(:);

% Filter BVAR monthly origins to quarter-ends
% dateNumeric_est was saved as the actual origin dates — use its length directly
nOrig_bv      = numel(bd.dateNumeric_est);
orig_dates_bv = datetime(bd.dateNumeric_est, 'ConvertFrom', 'datenum')';
qtr_mask_bv   = ismember(month(orig_dates_bv), [3 6 9 12]);
qtr_idx_bv    = find(qtr_mask_bv);
nQOrig_bv     = numel(qtr_idx_bv);
originDT_all{3} = dateshift(orig_dates_bv(qtr_mask_bv)', 'end', 'month');

zinf_bv = NaN(nHor, nQOrig_bv);
crps_bv = NaN(nHor, nQOrig_bv);
for qi = 1:nQOrig_bv
    t_idx = qtr_idx_bv(qi);
    for ih = 1:nHor
        k         = cfg.eval_horizons(ih);
        h_monthly = 3*k + 1;
        if h_monthly > size(pred_q_bv, 3), continue; end
        q_vals = squeeze(pred_q_bv(t_idx, :, h_monthly));
        if all(isnan(q_vals)), continue; end
        actual_h = actual_bv(t_idx, h_monthly);
        if isnan(actual_h), continue; end
        pit = interp1(q_vals(:), quant_bv, actual_h, 'linear', 'extrap');
        zinf_bv(ih, qi) = min(max(pit, 0), 1);
        % CRPS: generate samples from smoothed quantile distribution
        samp_bv = QR_sm(q_vals(:)', quant_bv(:)');
        crps_bv(ih, qi) = crps(samp_bv(:), actual_h, 1);
    end
end
zinf_all{3} = zinf_bv;
crps_all{3} = crps_bv;

clearvars bd pred_q_bv actual_bv quant_bv zinf_bv crps_bv samp_bv ...
          ih qi k h_monthly q_vals actual_h pit qtr_mask_bv qtr_idx_bv ...
          orig_dates_bv nOrig_bv nQOrig_bv;

%% ── Apply common evaluation window ──────────────────────────────────────
if ~isempty(cfg.eval_start) && ~isempty(cfg.eval_end)
    for m = 1:nModels
        keep = originDT_all{m} >= cfg.eval_start & originDT_all{m} <= cfg.eval_end;
        zinf_all{m}     = zinf_all{m}(:, keep);
        crps_all{m}     = crps_all{m}(:, keep);
        originDT_all{m} = originDT_all{m}(keep);
    end
end

% No Covid exclusion for inflation

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
        if numel(z_row) < 5, continue; end   % skip if too few obs

        % KS
        [results.(mod).ksinf(i), results.(mod).ksinfpv(i)] = kstestu(z_row);

        % Rossi-Sekhposyan (2019)
        [cv_i, stat_i, logic_i, KS_vec(i), CVM_vec(i)] = rs_test(z_row, rvec_rs);
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

        % CRPS mean
        crps_row = crps_all{m}(i, :);
        crps_row = crps_row(~isnan(crps_row));
        if ~isempty(crps_row)
            results.(mod).crps_mean(i) = mean(crps_row);
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
    KS   = results.(mod).KS_vec(:)';    % 1×H row vector
    CVM  = results.(mod).CVM_vec(:)';   % 1×H row vector

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

    results.(mod).gmm_ks       = mean(QVrejvecs,       1);
    results.(mod).gmm_cvm      = mean(CVMrejvecs,       1);
    results.(mod).gmm_ks_bonf  = mean(QVrejvecs_bonf,  1);
    results.(mod).gmm_cvm_bonf = mean(CVMrejvecs_bonf, 1);

    fprintf('GMM done: %s\n', mod);
end

%% ════════════════════════════════════════════════════════════════════════
%%  SAVE COMPARISON RESULTS
%%
%%  One Excel file: each sheet = one test statistic,
%%  columns = fanchart | qr | bvar, rows = horizons.
%% ════════════════════════════════════════════════════════════════════════

horiz    = cfg.eval_horizons(:);
filename = fullfile(evalDir, 'comparison_fcst_eval_inflation.xlsx');

% Helper: build side-by-side table for a given field
buildTab = @(field) table(horiz, ...
    results.fanchart.(field), results.qr.(field), results.bvar.(field), ...
    'VariableNames', {'horizon', 'fanchart', 'qr', 'bvar'});

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

% GMM summary (one row per model)
gmm_tab = table(models', ...
    [results.fanchart.gmm_ks;  results.qr.gmm_ks;  results.bvar.gmm_ks ], ...
    [results.fanchart.gmm_cvm; results.qr.gmm_cvm; results.bvar.gmm_cvm], ...
    'VariableNames', {'model','GMM_KS_std_bonf_w_wi','GMM_CVM_std_bonf_w_wi'});
writetable(gmm_tab, filename, 'Sheet', 'GMM_summary');

%% ════════════════════════════════════════════════════════════════════════
%%  PRINT SUMMARY TABLES  (one per model)
%% ════════════════════════════════════════════════════════════════════════

model_labels = {'Fan Chart (BOE)', 'QR (RASS)', 'BVAR'};

% Horizons to display
sel_h   = [0 1 2 3 4 8 12];
sel_idx = arrayfun(@(h) find(cfg.eval_horizons == h, 1), sel_h);

% F/R converter: logic arrays are stored as [1%, 5%, 10%] columns
% (0 = fail to reject null → 'F',  1 = reject → 'R',  NaN → '-')
frcode = @(v) char(70*(~isnan(v) & v < 0.5) + ...
                   82*(~isnan(v) & v >= 0.5) + ...
                   45*isnan(v));

hdr = ['  h   |  KS pval  | Knu pval  |  MW pval  | CRPS mean' ...
       ' | KS 10% | KS  5% | KS  1% | CvM 10% | CvM  5% | CvM  1%'];
sep = repmat('-', 1, numel(hdr));

for m = 1:nModels
    mod = models{m};
    fprintf('\n%s\n', sep);
    fprintf('  %s\n', model_labels{m});
    fprintf('%s\n', sep);
    fprintf('%s\n', hdr);
    fprintf('%s\n', sep);
    for ii = 1:numel(sel_idx)
        i  = sel_idx(ii);
        h  = cfg.eval_horizons(i);
        ks  = results.(mod).rs_ks_logic(i, :);   % columns: [1%  5%  10%]
        cvm = results.(mod).rs_cvm_logic(i, :);
        fprintf('  %2d  |  %7.4f  |  %7.4f  |  %7.4f  | %9.4f |   %s    |   %s    |   %s    |   %s     |   %s     |   %s\n', ...
            h, ...
            results.(mod).ksinfpv(i), ...
            results.(mod).knueppel_pval(i), ...
            results.(mod).MW_pval(i), ...
            results.(mod).crps_mean(i), ...
            frcode(ks(3)),  frcode(ks(2)),  frcode(ks(1)), ...   % 10%, 5%, 1%
            frcode(cvm(3)), frcode(cvm(2)), frcode(cvm(1)));     % 10%, 5%, 1%
    end
    fprintf('%s\n', sep);
end

%% ════════════════════════════════════════════════════════════════════════
%%  CHART: all tests × all models
%% ════════════════════════════════════════════════════════════════════════

colors     = [0.12 0.47 0.71;   % blue  – Fan Chart
              0.84 0.15 0.16;   % red   – QR
              0.17 0.63 0.17];  % green – BVAR
linestyles = {'-', '--', ':'};
markers    = {'o', 's', '^'};
lw         = 1.8;
h_vec      = cfg.eval_horizons(:);

test_fields = {'ksinfpv',   'knueppel_pval',          'MW_pval',                    'crps_mean'};
test_titles = {'KS p-value','Knueppel (2015) p-value','Mitchell-Weale (2023) p-value','CRPS mean'};
is_pval     = logical([1 1 1 0]);

fig = figure('Name','Forecast Evaluation — CPI Inflation','NumberTitle','off', ...
             'Position',[80 80 1200 820]);

for t = 1:4
    ax = subplot(2, 2, t);
    hold(ax, 'on');

    for m = 1:nModels
        mod = models{m};
        y   = results.(mod).(test_fields{t});
        plot(ax, h_vec, y, ...
             'Color',     colors(m,:), ...
             'LineStyle', linestyles{m}, ...
             'LineWidth', lw, ...
             'Marker',    markers{m}, ...
             'MarkerSize', 5, ...
             'MarkerFaceColor', colors(m,:), ...
             'DisplayName', model_labels{m});
    end

    % significance lines for p-value panels
    if is_pval(t)
        yline(ax, 0.10, 'Color',[0.4 0.4 0.4], 'LineStyle',':', ...
              'LineWidth',1.0, 'Label','10%', 'LabelVerticalAlignment','bottom', ...
              'HandleVisibility','off');
        yline(ax, 0.05, 'Color',[0.2 0.2 0.2], 'LineStyle','--', ...
              'LineWidth',1.0, 'Label','5%',  'LabelVerticalAlignment','bottom', ...
              'HandleVisibility','off');
        ylim(ax, [0 1]);
    end

    xlabel(ax, 'Horizon (quarters ahead)', 'FontSize', 9);
    title(ax,  test_titles{t}, 'FontSize', 10, 'FontWeight','bold');
    set(ax, 'XTick', h_vec, 'Box','on', 'FontSize', 9);
    grid(ax, 'on');

    if t == 1   % single legend in top-left panel
        lgd = legend(ax, 'Location','best', 'FontSize', 9);
        lgd.Box = 'on';
    end

    hold(ax, 'off');
end

sgtitle('CPI Inflation — Forecast Evaluation', 'FontSize', 13, 'FontWeight','bold');

% Save
chartFile = fullfile(evalDir, 'comparison_fcst_eval_inflation.png');
exportgraphics(fig, chartFile, 'Resolution', 300);
fprintf('\nChart saved to %s\n', chartFile);

fprintf('\n── DONE: comparison saved to %s ──\n', filename);
