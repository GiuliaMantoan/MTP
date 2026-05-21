%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FORECAST EVALUATION — GDP GROWTH  (multi-model comparison)
%%
%%  Authors: David Aikman, Rhys Bidder, Simon Lloyd, Giulia Mantoan,
%%           Simone Maso, Aditya Mori, Matthew Tong
%%
%%  Tests: KS · Rossi-Sekhposyan (2019) · Berkowitz (2001) ·
%%         Knueppel (2015) · Mitchell-Weale (2023) · Galvao-Mantoan-Mitchell
%%
%%  Compares: Fan Chart · QR · BVAR
%%  Results are saved side-by-side for each test statistic.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% ── Paths (defined first so scriptDir is available throughout) ───────────
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
% Origins outside this window are treated as NaN in the PITs.
cfg.eval_start = datetime('30-Sep-2007', 'InputFormat', 'dd-MMM-yyyy');
cfg.eval_end   = datetime('31-Mar-2024', 'InputFormat', 'dd-MMM-yyyy');

% Covid exclusion window (same treatment as GDP QR and BVAR models)
% Origins whose date falls in [covidStart, covidEnd] are set to NaN in all models.
cfg.covidStart = datetime(2020, 3, 31);   % Q1 2020 (first excluded quarter)
cfg.covidEnd   = datetime(2022, 3, 31);   % Q1 2022 (annualised rate still spans Covid)

% Fan chart settings
cfg.fanchart.start_date = datetime('30-Sep-2007', 'InputFormat', 'dd-MMM-yyyy');
cfg.fanchart.end_date   = datetime('31-Mar-2024', 'InputFormat', 'dd-MMM-yyyy');
cfg.fanchart.end_fcst   = datetime('30-Mar-2027', 'InputFormat', 'dd-MMM-yyyy');
cfg.fanchart.covid_date = datetime('30-Jun-2020', 'InputFormat', 'dd-MMM-yyyy');

% Root directory that boefsctdata_growth uses via cd() to locate
% Data\kcl_data\gdp_growth_projection_parameters_mpc.xlsx
% Set this to the folder that contains the Data\ subfolder.
cfg.paths.boe_root     = scriptDir;   % ← adjust if the BoE Excel is elsewhere
cfg.paths.fanchart_raw = fullfile(scriptDir, 'GaRDataRaw_quarterly_BIS_march.xlsx');
cfg.paths.qr_pred_q    = fullfile(outDir, 'predicted_quantiles_gdp_OOS.mat');
cfg.paths.qr_actual    = fullfile(outDir, 'actual_gdp_yoy_OOS.mat');
cfg.paths.bvar         = fullfile(outDir, 'BVAR', 'BVAR_gdp_pred_q.mat');

%% ── Global settings ──────────────────────────────────────────────────────
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth', 1);
rng(0, 'twister');   % specify generator explicitly so legacy mode doesn't block this

addpath(fullfile(scriptDir, 'intermediate_codes'));
addpath(fullfile(scriptDir, 'functions'));
addpath(fullfile(scriptDir, 'functions', 'azzalini'));

nHor    = numel(cfg.eval_horizons);
models  = {'fanchart', 'qr', 'bvar'};
nModels = numel(models);

%% ════════════════════════════════════════════════════════════════════════
%%  LOAD ALL MODELS AND COMPUTE PITs
%%
%%  zgrowth_all{m} — nHor × nOrigins_m  matrix of PITs for model m.
%%  Row ih = cfg.eval_horizons(ih) quarters ahead.
%%  Columns not in cfg.eval_start/end window are left as NaN.
%% ════════════════════════════════════════════════════════════════════════

zgrowth_all   = cell(nModels, 1);
originDT_all  = cell(nModels, 1);   % datetime of each origin (for window filter)

%% ── (1) Fan chart ────────────────────────────────────────────────────────
fprintf('Loading fan chart...\n');

% Variables required by boefsctdata_growth and actualdata
fullFileName = cfg.paths.fanchart_raw;   % used by actualdata
varnames     = {'g4rgdp'};
ctrynames    = {'UK'};
momentlist   = {'fcstmean', 'fcststdev', 'fcstskew'};
outputFolder = evalDir;

start_date = cfg.fanchart.start_date;
end_date   = cfg.fanchart.end_date;
end_fcst   = cfg.fanchart.end_fcst;     % required by boefsctdata_growth
covid_date = cfg.fanchart.covid_date;

% boefsctdata_growth ends with clearvars -except ..., wiping most of the
% workspace including ws_backup itself. We save using outputFolder directly
% (it IS in the exception list and survives), then rebuild the path from it.
save(fullfile(outputFolder, 'ws_backup_temp.mat'));

cd(cfg.paths.boe_root);
boefsctdata_growth;   % loads mtestdata; then wipes workspace via clearvars
% ← only mtestdata + {start_date end_date covid_date ctrynames varnames
%   fullFileName momentlist outputFolder covid_end_date} survive here

mtestdata_fc = mtestdata;                           % grab before restore
load(fullfile(outputFolder, 'ws_backup_temp.mat')); % outputFolder still alive
delete(fullfile(evalDir, 'ws_backup_temp.mat'));    % evalDir restored by load
mtestdata = mtestdata_fc;                           % put fanchart data back
clear mtestdata_fc;

% actualdata has no clearvars — safe to call directly.
% It needs: fullFileName, varnames, ctrynames, start_date, end_date.
actualdata;           % loads actualvar  (13 × nOrigins)

% Covid column drop
covid_col = (year(covid_date) - year(start_date)) * 4 + ...
            (ceil(month(covid_date)/3) - ceil(month(start_date)/3));
if size(actualvar, 2) >= (covid_col + 1)
    actualvar(:, covid_col + 1) = [];
end
if size(mtestdata, 2) >= (covid_col * 3 + 1)
    mtestdata(:, (covid_col*3+1):(covid_col*3+3)) = [];
end

modevar = mtestdata(:, 1:3:end);
meanvar = mtestdata(:, 2:3:end);
dispvar = mtestdata(:, 3:3:end);

% Build quarter dates and drop the same Q2-2020 column removed from actualvar
quarterDates_full = start_date : calmonths(3) : end_date;
quarterDates = [quarterDates_full(1:covid_col), quarterDates_full(covid_col+2:end)];
nOrig_fc     = size(meanvar, 2);
quarterDates = quarterDates(1 : nOrig_fc);   % guard against any length mismatch
originDT_all{1} = quarterDates(:);

zgrowth_fc = NaN(size(actualvar));
for i = 1:size(meanvar, 1)
    for j = 1:nOrig_fc
        if meanvar(i,j) - modevar(i,j) ~= 0
            gam = stdtogam(meanvar(i,j), modevar(i,j), dispvar(i,j));
            sig = mom2g(dispvar(i,j), gam);
            [m, s1, s2] = momtopar(meanvar(i,j), modevar(i,j), sig);
            if ~isnan(actualvar(i,j))
                zgrowth_fc(i,j) = integral(@(y) ftp(y,m,s1,s2), -3, actualvar(i,j));
            end
        else
            if ~isnan(actualvar(i,j))
                zgrowth_fc(i,j) = normcdf((actualvar(i,j) - meanvar(i,j)) / dispvar(i,j));
            end
        end
    end
end
zgrowth_fc = zgrowth_fc(:, ~all(isnan(zgrowth_fc), 1));

% Restrict to eval horizons (rows = horizons in fan chart are 0..H-1)
if size(zgrowth_fc, 1) >= nHor
    zgrowth_all{1} = zgrowth_fc(cfg.eval_horizons + 1, :);
else
    zgrowth_all{1} = zgrowth_fc;
end

clearvars modevar meanvar dispvar zgrowth_fc actualvar mtestdata ...
          gam sig i j covid_col nOrig_fc;

%% ── (2) QR ───────────────────────────────────────────────────────────────
fprintf('Loading QR model...\n');

qd = load(cfg.paths.qr_pred_q, 'pred_q');
ad = load(cfg.paths.qr_actual,  'actual_var', 'idx_est', 'dateNumeric_full');
pred_q_qr  = qd.pred_q;       % nOrigins × Qn × horizons
actual_qr  = ad.actual_var;
idx_est_qr = ad.idx_est;
dates_qr   = ad.dateNumeric_full;
quant_qr   = (0.05:0.05:0.95)';

nOrig_qr        = size(pred_q_qr, 1);
originDT_all{2} = datetime(dates_qr(idx_est_qr : idx_est_qr + nOrig_qr - 1), ...
                            'ConvertFrom','datenum')';

zgrowth_qr = NaN(nHor, nOrig_qr);
for ih = 1:nHor
    h = cfg.eval_horizons(ih) + 1;          % pred_q index (h=1 = nowcast)
    if h > size(pred_q_qr, 3), continue; end
    for t = 1:nOrig_qr
        q_vals = squeeze(pred_q_qr(t, :, h));
        if all(isnan(q_vals)), continue; end
        tgt = (idx_est_qr + t - 1) + (h - 1);
        if tgt > length(actual_qr), continue; end
        actual_h = actual_qr(tgt);
        if isnan(actual_h), continue; end
        pit = interp1(q_vals(:), quant_qr, actual_h, 'linear', 'extrap');
        zgrowth_qr(ih, t) = min(max(pit, 0), 1);
    end
end
zgrowth_all{2} = zgrowth_qr;
clearvars qd ad pred_q_qr actual_qr idx_est_qr dates_qr quant_qr zgrowth_qr ...
          ih h t q_vals tgt actual_h pit nOrig_qr;

%% ── (3) BVAR ─────────────────────────────────────────────────────────────
fprintf('Loading BVAR...\n');

bd = load(cfg.paths.bvar, 'pred_q', 'actualvar', 'dateNumeric_est', 'cfg');
pred_q_bv  = bd.pred_q;       % nOrigins × Qn × horizons
actual_bv  = bd.actualvar;    % nOrigins × horizons  (already aligned)
quant_bv   = bd.cfg.quantiles(:);

nOrig_bv        = size(pred_q_bv, 1);
originDT_all{3} = datetime(bd.dateNumeric_est, 'ConvertFrom','datenum')';

zgrowth_bv = NaN(nHor, nOrig_bv);
for ih = 1:nHor
    h = cfg.eval_horizons(ih) + 1;
    if h > size(pred_q_bv, 3), continue; end
    for t = 1:nOrig_bv
        q_vals = squeeze(pred_q_bv(t, :, h));
        if all(isnan(q_vals)), continue; end
        actual_h = actual_bv(t, h);
        if isnan(actual_h), continue; end
        pit = interp1(q_vals(:), quant_bv, actual_h, 'linear', 'extrap');
        zgrowth_bv(ih, t) = min(max(pit, 0), 1);
    end
end
zgrowth_all{3} = zgrowth_bv;
clearvars bd pred_q_bv actual_bv quant_bv zgrowth_bv ...
          ih h t q_vals actual_h pit nOrig_bv;

%% ── Apply common evaluation window (if set) ──────────────────────────────
if ~isempty(cfg.eval_start) && ~isempty(cfg.eval_end)
    for m = 1:nModels
        keep = originDT_all{m} >= cfg.eval_start & originDT_all{m} <= cfg.eval_end;
        zgrowth_all{m}  = zgrowth_all{m}(:, keep);
        originDT_all{m} = originDT_all{m}(keep);
    end
end

%% ── Apply Covid exclusion window (same treatment as GDP QR and BVAR) ─────
% Origins in [covidStart, covidEnd] are set to NaN in all three models.
% For QR and BVAR this is already NaN from the models themselves;
% this step enforces the same exclusion uniformly on the fan chart too.
for m = 1:nModels
    covid_mask = originDT_all{m} >= cfg.covidStart & originDT_all{m} <= cfg.covidEnd;
    zgrowth_all{m}(:, covid_mask) = NaN;
end

%% ════════════════════════════════════════════════════════════════════════
%%  RUN TESTS FOR EACH MODEL
%% ════════════════════════════════════════════════════════════════════════

% Pre-allocate result structs (one entry per model)
results = struct();
for m = 1:nModels
    results.(models{m}).ksgrowth     = NaN(nHor, 1);
    results.(models{m}).ksgrowthpv   = NaN(nHor, 1);
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
end

lags     = -1;    % automatic lag selection (Knueppel / MW)
prewhite =  0;
z_u = 0.95;  z_l = 0.05;  df = 5;   % MW parameters
rvec_rs = linspace(0, 1, 1000);

for m = 1:nModels
    mod = models{m};
    zgrowth = zgrowth_all{m};
    fprintf('\nRunning tests: %s ...\n', mod);

    KS_vec  = zeros(nHor, 1);
    CVM_vec = zeros(nHor, 1);

    for i = 1:nHor
        z_row = zgrowth(i, :)';
        z_row = z_row(~isnan(z_row));
        if numel(z_row) < 5, continue; end   % skip if too few obs

        % KS
        [results.(mod).ksgrowth(i), results.(mod).ksgrowthpv(i)] = kstestu(z_row);

        % Rossi-Sekhposyan
        [cv_i, stat_i, logic_i, KS_vec(i), CVM_vec(i)] = rs_test(z_row, rvec_rs);
        results.(mod).rs_ks_test(i)   = stat_i.ks_stat;
        results.(mod).rs_cvm_test(i)  = stat_i.cvm_stat;
        results.(mod).rs_ks_logic(i,:)  = logic_i.ks_array;
        results.(mod).rs_cvm_logic(i,:) = logic_i.cvm_array;

        % Berkowitz
        [results.(mod).bert1(i), results.(mod).bert2(i), ...
         results.(mod).berK1(i), results.(mod).berK2(i)] = berk(zgrowth(i, :));

        % Knueppel
        [results.(mod).knueppel_stat(i), results.(mod).knueppel_pval(i)] = ...
            alpha0_1234_NW(zgrowth(i,:), lags, prewhite);

        % Mitchell-Weale
        stat_u = MW_alpha0_1234_NW(zgrowth(i,:), lags, prewhite, z_l, z_u);
        stat_f = freqTestInCensoredRegion(zgrowth(i,:), lags, prewhite, z_l, z_u);
        results.(mod).MW_stat(i) = stat_u + stat_f;
        results.(mod).MW_pval(i) = 1 - chi2cdf(results.(mod).MW_stat(i), df);
    end

    results.(mod).KS_vec  = KS_vec;
    results.(mod).CVM_vec = CVM_vec;
end

%% Galvao-Mantoan-Mitchell  (run once per model — computationally intensive)
MC     = 1000;
bootMC = 1000;
rng(bootMC, 'twister');   % modern equivalent of rand/randn('seed', bootMC)
rvec_gmm = 0:0.001:1;

for m = 1:nModels
    mod     = models{m};
    zgrowth = zgrowth_all{m};
    z       = zgrowth';
    P       = size(z, 1);
    el      = floor(P^(1/4));
    Hz      = size(z, 2);
    KS      = results.(mod).KS_vec(:)';   % size_statistic_h2 expects 1×H row vector
    CVM     = results.(mod).CVM_vec(:)';  % size_statistic_h2 expects 1×H row vector

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

    results.(mod).gmm_ks  = mean(QVrejvecs,       1);
    results.(mod).gmm_cvm = mean(CVMrejvecs,       1);
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
filename = fullfile(evalDir, 'comparison_fcst_eval_growth.xlsx');

% Helper: build side-by-side table for a given field
buildTab = @(field) table(horiz, ...
    results.fanchart.(field), results.qr.(field), results.bvar.(field), ...
    'VariableNames', {'horizon', 'fanchart', 'qr', 'bvar'});

writetable(buildTab('ksgrowth'),      filename, 'Sheet', 'KS_stat');
writetable(buildTab('ksgrowthpv'),    filename, 'Sheet', 'KS_pval');
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

% GMM summary (one row per model)
gmm_tab = table(models', ...
    [results.fanchart.gmm_ks;  results.qr.gmm_ks;  results.bvar.gmm_ks ], ...
    [results.fanchart.gmm_cvm; results.qr.gmm_cvm; results.bvar.gmm_cvm], ...
    'VariableNames', {'model','GMM_KS_std_bonf_w_wi','GMM_CVM_std_bonf_w_wi'});
writetable(gmm_tab, filename, 'Sheet', 'GMM_summary');

fprintf('\n── DONE: comparison saved to %s ──\n', filename);
