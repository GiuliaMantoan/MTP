%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  INFLATION AT-RISK  —  Recursive OOS Quantile LP + Model/Spec Selection
%%  (Monthly, Multi-Specification)
%
%  Authors : Aikman, Bidder, Lloyd, Mantoan, Maso, Mori, Tong
%  Updated : April 2026
%
%  PIPELINE
%    1. Generate all variable-category combinations  (combo_specifications)
%    2. For each spec: load data, estimate QR, build predicted quantiles
%    3. Select best spec via Weighted Interval Score (WIS)
%    4. For selected spec only: fit distribution (skew-t / semi-param / TPN)
%    5. Save .mat outputs; export moments to Excel
%    6. Figures: rolling fan / forward fan / HD / López-Salido (1m,1y,2y)
%    7. Export last-origin skew-t parameters for sharing
%
%  MODEL SELECTION   (cfg.model_selection)
%    2 = QR + Skew-t (Azzalini-Capitanio)
%    3 = QR + Semi-parametric (Mitchell-Poon-Zhu)
%    4 = QR + Two-piece Normal
%
%  SPEC SELECTION    (cfg.use_best_spec)
%    1 = auto: pick spec with lowest average WIS across all horizons
%    0 = manual: use cfg.specplot
%
%  Add alternatives to each cfg.var.* cell to search over more specs, e.g.:
%    cfg.var.slack = {'delta1vu_tightness', 'ugap_hp_filter_lambda_129600'};
%
%  NOTE: Set cfg.bst.nboot = 5000 for production runs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% ── PATHS ────────────────────────────────────────────────────────────────
scriptDir = fileparts(mfilename('fullpath'));

dataFile = fullfile(scriptDir, 'IaRDataRaw_monthly_M.xlsx');
outDir   = fullfile(scriptDir, 'Outputs');
figDir   = fullfile(outDir, 'econ_interpretation_charts');
fanDir   = fullfile(outDir, 'predictive_densities');
sktDir   = fullfile(outDir, 'sktparam');
wisDir   = fullfile(outDir, 'wis');

mkdirs({outDir, figDir, fanDir, sktDir, wisDir});

addpath(fullfile(scriptDir, 'intermediate_codes'));
addpath(fullfile(scriptDir, 'functions'));
addpath(fullfile(scriptDir, 'functions', 'azzalini'));
addpath(fullfile(scriptDir, 'functions', 'CRPS'));
addpath(fullfile(scriptDir, 'functions', 'Simon_qreg'));
addpath(fullfile(scriptDir, 'functions', 'two_piece_normal'));

%% ── CONFIGURATION ────────────────────────────────────────────────────────

% Sample
cfg.startT   = datenum(1990, 2, 1);
cfg.endT     = datenum(2026, 3, 1);
cfg.startEst = datenum(2004, 1, 1);    % first forecast origin

% Dependent variable (single element; changing it invalidates WIS comparisons)
cfg.var.dep = {'g4cpi'};

% Variable categories — list alternatives in each cell to search over specs.
% The ndgrid below generates every combination automatically.
% e.g. two options in slack × two in supply = 4 specs total.

cfg.var.persist = {'avgcpi'};                     % inflation persistence
%                  'avg12Infl'                    % (alternative: 12m avg CPI)

cfg.var.expect  = {'infl1_m'};                    % inflation expectations
%                  'infl2_m'                       % (alternative: 2y ahead expectations)

cfg.var.slack   = {'delta1vu_tightness'};         % economic slack
%                  'ugap_hp_filter_lambda_129600'  % (alternative: unemployment gap HP)
%                  'ugap_kalman_filter'            % (alternative: unemployment gap Kalman)
%                  'vu_tightness'                  % (alternative: v/u ratio level)

cfg.var.supply  = {'yoy_growth_import_deflator'}; % supply / external
%                  'g4oil'                         % (alternative: oil price inflation)
%                  'g4cpicore_global'              % (alternative: G7 core inflation)

cfg.var.fci     = {'bond_spread'};                % financial conditions
%                  'market_vol_uk'                 % (alternative: market volatility)

% Distribution model
%   2 = Skew-t  |  3 = Semi-parametric  |  4 = Two-piece Normal
cfg.model_selection = 2;
cfg.modellist       = {'ols', 'skewt', 'semi-param', 'two-piece-normal'};

% Spec selection
cfg.use_best_spec = 1;       % 1 = auto (lowest avg WIS);  0 = manual
cfg.specplot      = 1;       % used when use_best_spec = 0
cfg.horizons_wis  = [1 13 25]; % horizons for the per-horizon WIS table

% Covid dummy (off for inflation)
cfg.covid = 0;

% Model
cfg.horizons  = 37;            % h=1 current; h=2..37 → 1..36 months ahead
cfg.quantiles = 0.05:0.05:0.95;

% Bootstrap  *** SET nboot = 5000 FOR PRODUCTION RUNS ***
cfg.bst.nboot     = 10;
cfg.bst.blocksize = 24;
cfg.bst.ci        = 68;

% Plotting
cfg.hPlot   = [2, 4, 13];             % horizons for rolling fan + HD charts
cfg.qDecomp = [0.25, 0.50, 0.75];     % quantiles for HD charts
cfg.hor_ls  = [2, 13, 25];            % López-Salido horizons (1m, 1y, 2y)

%% ── GLOBAL PLOT DEFAULTS ─────────────────────────────────────────────────
set(0, 'defaultAxesFontName',      'Times');
set(0, 'defaultAxesLineStyleOrder', '-|--|:');
set(0, 'defaultLineLineWidth',       1);
rng('default');

%% ════════════════════════════════════════════════════════════════════════
%%  1.  GENERATE SPEC COMBINATIONS
%% ════════════════════════════════════════════════════════════════════════

[i0,i1,i2,i3,i4,i5] = ndgrid( ...
    1:numel(cfg.var.dep),     1:numel(cfg.var.persist), ...
    1:numel(cfg.var.expect),  1:numel(cfg.var.slack),   ...
    1:numel(cfg.var.supply),  1:numel(cfg.var.fci));

vars = { cfg.var.dep(i0(:)),     cfg.var.persist(i1(:)), ...
         cfg.var.expect(i2(:)),  cfg.var.slack(i3(:)),   ...
         cfg.var.supply(i4(:)),  cfg.var.fci(i5(:)) };
for k = 1:numel(vars)
    if size(vars{k}, 2) > 1, vars{k} = vars{k}.'; end
end

combo_specifications = [vars{:}];          % nSpec × (1 dep + 5 predictors)
nSpec = size(combo_specifications, 1);
nPred = size(combo_specifications, 2) - 1; % predictors (excl. dep var)
V     = nPred + 1;                         % const + predictors

spec_names = arrayfun(@(s) strjoin(combo_specifications(s,:), ' | '), ...
             (1:nSpec)', 'uni', false);
fprintf('Total specifications: %d\n', nSpec);

%% ════════════════════════════════════════════════════════════════════════
%%  2.  PRE-ALLOCATE ARRAYS
%% ════════════════════════════════════════════════════════════════════════

Qn   = numel(cfg.quantiles);
dv1  = datevec(cfg.startEst); dv2 = datevec(cfg.endT);
nOrigins = (dv2(1)-dv1(1))*12 + (dv2(2)-dv1(2)) + 1;

pred_q      = NaN(nOrigins, Qn, cfg.horizons, nSpec);
coeffqr     = NaN(V, Qn, cfg.horizons, nOrigins, nSpec);
bootstrapqr = NaN(V, Qn, cfg.horizons, cfg.bst.nboot, nSpec);
avg_wis     = NaN(nSpec, 1);

% Shared across specs (set once during first iteration)
idx_est          = [];
last_origin      = [];
dateNumeric_full = [];
actualvar        = [];   % nOrigins × horizons  (actual dep-var for WIS)

%% ════════════════════════════════════════════════════════════════════════
%%  3.  SPEC LOOP  —  QR ESTIMATION + WIS
%% ════════════════════════════════════════════════════════════════════════

tic;
for spec = 1:nSpec

    fprintf('\n── Spec %d/%d: %s\n', spec, nSpec, spec_names{spec});

    %% ── Load data ────────────────────────────────────────────────────────
    varnames     = combo_specifications(spec,:);
    lagctryvar   = ones(1, numel(varnames));
    lagglobalvar = [];
    lags         = lagctryvar;
    startT       = cfg.startT;      endT         = cfg.endT;
    ctrynames    = {'UK'};          onlyuk       = 1;
    fullFileName = dataFile;        covid        = cfg.covid;
    dummyvarname = {'covid'};       quantilelevels = cfg.quantiles;
    horizons     = cfg.horizons;    bstOptions   = cfg.bst;

    dataloading_monthly;

    %% ── First-iteration setup (spec-independent quantities) ─────────────
    if spec == 1
        idx_est          = find(dateNumeric == cfg.startEst, 1);
        last_origin      = numel(dateNumeric);
        dateNumeric_full = dateNumeric;
        assert(~isempty(idx_est), 'cfg.startEst not found in data.');

        % actualvar(t,h) = dep-var value at the h-step-ahead target of origin t
        act_long = countryData.UK.(cfg.var.dep{1});
        act_filt = act_long(idx_est:end);
        n_act    = numel(act_filt);
        interm   = NaN(n_act, n_act);
        for ii = 1:n_act
            interm(1:(n_act-ii+1), ii) = act_filt(ii:end);
        end
        actualvar = interm(1:cfg.horizons, :)';   % nOrigins × horizons
    end

    %% ── Extract predictors ───────────────────────────────────────────────
    vlag    = strcat('l1', varnames);
    depvar  = countryData.UK.(vlag{1});
    explvar = table2array(countryData.UK(:, vlag(2:end)));

    %% ── Recursive QR loop ───────────────────────────────────────────────
    for endtime = idx_est : last_origin

        t = endtime - idx_est + 1;
        X = explvar(1:endtime, :);
        y = depvar( 1:endtime, :);

        Y_LP = NaN(size(X,1), 1, cfg.horizons);
        for h = 1:cfg.horizons
            if size(X,1) > h
                Y_LP(1:size(X,1)-h, 1, h) = y(1+h:end);
            end
        end

        doBootstrap = (endtime == last_origin);
        [bQR, bQRbst] = qfe_qr_local_projection_SL_final( ...
            Y_LP, [ones(size(X,1),1), X], cfg.quantiles, (1:cfg.horizons)', ...
            doBootstrap, cfg.bst, 0, 0);

        coeffqr(:,:,:,t,spec) = bQR;
        if doBootstrap, bootstrapqr(:,:,:,:,spec) = bQRbst; end

        x_now = [1, X(endtime,:)];
        for h = 1:cfg.horizons
            pred_q(t,:,h,spec) = x_now * bQR(:,:,h);
        end
    end

    pred_q(:,:,:,spec) = sort(pred_q(:,:,:,spec), 2);   % enforce monotonicity

    %% ── WIS (vectorised over quantiles) ─────────────────────────────────
    ql3   = reshape(cfg.quantiles, 1, Qn, 1);
    av3   = repmat(reshape(actualvar, nOrigins, 1, cfg.horizons), [1, Qn, 1]);
    u3    = av3 - squeeze(pred_q(:,:,:,spec));
    check = max(ql3 .* u3, (ql3-1) .* u3);
    avg_wis(spec) = mean(check(~isnan(av3(:))), 'omitnan');
    fprintf('  WIS = %.4f\n', avg_wis(spec));

end
fprintf('\nEstimation wall-clock: %.1f s\n', toc);

%% ════════════════════════════════════════════════════════════════════════
%%  4.  SPEC SELECTION
%% ════════════════════════════════════════════════════════════════════════

[~, idx_min] = min(avg_wis);
fprintf('\nBest spec (avg WIS = %.4f): %s\n', avg_wis(idx_min), spec_names{idx_min});

% Per-horizon WIS table
wis_by_hor = NaN(nSpec, cfg.horizons);
for s = 1:nSpec
    for h = 1:cfg.horizons
        pq_h = squeeze(pred_q(:,:,h,s));   % nOrigins × Qn
        av_h = actualvar(:,h);
        mask = ~isnan(av_h);
        if any(mask)
            u_h  = repmat(av_h(mask), [1,Qn]) - pq_h(mask,:);
            ql_  = repmat(cfg.quantiles,    [sum(mask),1]);
            wis_by_hor(s,h) = mean(max(ql_.*u_h, (ql_-1).*u_h), 'all');
        end
    end
end
[~, idx_min_hor] = min(wis_by_hor(:, cfg.horizons_wis), [], 1);

% Save WIS tables to Excel
xlsWIS = fullfile(wisDir, sprintf('wis_%s_inflation.xlsx', ...
    cfg.modellist{cfg.model_selection}));
writetable(table(spec_names, avg_wis, 'VariableNames',{'Specification','AvgWIS'}), ...
    xlsWIS, 'Sheet','avg_wis', 'WriteMode','overwrite');
writetable(table(cfg.horizons_wis(:), spec_names(idx_min_hor(:)), ...
    'VariableNames',{'Horizon','BestSpec'}), ...
    xlsWIS, 'Sheet','best_by_horizon', 'WriteMode','append');

% Determine which spec to use for all downstream analysis
if cfg.use_best_spec
    spec_to_use = idx_min;
else
    spec_to_use = cfg.specplot;
end
fprintf('Using spec %d: %s\n', spec_to_use, spec_names{spec_to_use});

%% ── Reduce arrays to spec_to_use ────────────────────────────────────────
pred_q      = squeeze(pred_q(:,:,:,spec_to_use));         % nO × Qn × H
coeffqr     = squeeze(coeffqr(:,:,:,:,spec_to_use));      % V  × Qn × H × nO
bootstrapqr = squeeze(bootstrapqr(:,:,:,:,spec_to_use));  % V  × Qn × H × nboot

%% ── Re-load data for spec_to_use ────────────────────────────────────────
varnames     = combo_specifications(spec_to_use,:);
lagctryvar   = ones(1, numel(varnames));
lagglobalvar = [];   lags = lagctryvar;
startT       = cfg.startT;      endT         = cfg.endT;
ctrynames    = {'UK'};          onlyuk       = 1;
fullFileName = dataFile;        covid        = cfg.covid;
dummyvarname = {'covid'};       quantilelevels = cfg.quantiles;
horizons     = cfg.horizons;    bstOptions   = cfg.bst;
dataloading_monthly;

vlag    = strcat('l1', varnames);
depvar  = countryData.UK.(vlag{1});
explvar = table2array(countryData.UK(:, vlag(2:end)));
actual_var = countryData.UK.(cfg.var.dep{1});

StartEstDT    = datetime(cfg.startEst,              'ConvertFrom','datenum');
lastOriDT     = datetime(dateNumeric_full(last_origin), 'ConvertFrom','datenum');
months_origin = StartEstDT : calmonths(1) : lastOriDT;   % 1 × nOrigins

Xorig   = explvar(idx_est:last_origin, :);
Xfull   = [ones(nOrigins,1), Xorig];   % nOrigins × V
vLabels = [{'Constant'}, ...
    cellfun(@(s) strrep(s,'_',' '), combo_specifications(spec_to_use,2:end), 'uni',0)];

save(fullfile(outDir,'explanatoryvar_inflation_OOS.mat'), 'explvar');

%% ════════════════════════════════════════════════════════════════════════
%%  5.  DISTRIBUTION FITTING  (for spec_to_use only)
%% ════════════════════════════════════════════════════════════════════════

switch cfg.model_selection

    case 2  %────────────────── Skew-t (Azzalini-Capitanio) ───────────────
        lc_skt = NaN(nOrigins, cfg.horizons);   sc_skt = NaN(nOrigins, cfg.horizons);
        sh_skt = NaN(nOrigins, cfg.horizons);   df_skt = NaN(nOrigins, cfg.horizons);
        fcstmean  = NaN(nOrigins, cfg.horizons);
        fcststdev = NaN(nOrigins, cfg.horizons);
        fcstskew  = NaN(nOrigins, cfg.horizons);

        for t = 1:nOrigins
            for h = 1:cfg.horizons
                [lc_skt(t,h), sc_skt(t,h), sh_skt(t,h), df_skt(t,h)] = ...
                    QuantilesInterpolation(pred_q(t,:,h), cfg.quantiles);
                [fcstmean(t,h), fcststdev(t,h), fcstskew(t,h), ~] = ...
                    moments_skewt_distr_sm_updated( ...
                        lc_skt(t,h), sc_skt(t,h), sh_skt(t,h), df_skt(t,h));
            end
        end

    case 3  %────────────────── Semi-parametric (Mitchell-Poon-Zhu) ───────
        semi_param_distr = NaN(nOrigins, 20000, cfg.horizons);
        empirical_cdf    = NaN(20001, 2, nOrigins, cfg.horizons);
        fcstmean  = NaN(nOrigins, cfg.horizons);
        fcststdev = NaN(nOrigins, cfg.horizons);
        fcstskew  = NaN(nOrigins, cfg.horizons);

        for h = 1:cfg.horizons
            semi_param_distr(:,:,h) = QR_sm(pred_q(:,:,h), cfg.quantiles);
            fcstmean(:,h)  = mean(semi_param_distr(:,:,h), 2);
            fcststdev(:,h) = std( semi_param_distr(:,:,h), 0, 2);
            fcstskew(:,h)  = skewness(semi_param_distr(:,:,h), 0, 2);
            for t = 1:nOrigins
                [f,x] = ecdf(semi_param_distr(t,:,h));
                blk = nan(20001,2);  blk(1:numel(x),:) = [x, f];
                empirical_cdf(:,:,t,h) = blk;
            end
        end

        % CRPS for selected spec
        crps_results = NaN(nOrigins, cfg.horizons);
        for h = 1:cfg.horizons
            for t = 1:nOrigins
                crps_results(t,h) = crps(semi_param_distr(t,:,h), actualvar(t,h), 2);
            end
        end
        fprintf('Avg CRPS (spec_to_use): %.4f\n', mean(crps_results(:),'omitnan'));

    case 4  %────────────────── Two-piece Normal ──────────────────────────
        param_tpn = NaN(nOrigins, 3, cfg.horizons);
        fcstmean  = NaN(nOrigins, cfg.horizons);
        fcststdev = NaN(nOrigins, cfg.horizons);
        fcstskew  = NaN(nOrigins, cfg.horizons);

        for t = 1:nOrigins
            fprintf('TPN fit: origin %d/%d\n', t, nOrigins);
            for h = 1:cfg.horizons
                param_tpn(t,:,h) = fit_tpn_to_quantiles_SM(pred_q(t,:,h), cfg.quantiles);
                fcstmean(t,h)    = two_part_normal_mean(    param_tpn(t,:,h));
                fcststdev(t,h)   = sqrt(two_part_normal_variance(param_tpn(t,:,h)));
                fcstskew(t,h)    = two_part_normal_skewness(param_tpn(t,:,h));
            end
        end

end

%% ════════════════════════════════════════════════════════════════════════
%%  6.  SAVE CORE OUTPUTS
%% ════════════════════════════════════════════════════════════════════════

save(fullfile(outDir,'actual_inflation_mom_OOS.mat'), ...
    'actual_var','dateNumeric_full','idx_est','last_origin','startT','endT');
save(fullfile(outDir,'qreg_results_inflation_OOS.mat'),        'coeffqr');
save(fullfile(outDir,'predicted_quantiles_inflation_OOS.mat'), 'pred_q');
save(fullfile(outDir,'bootstrap_results_inflation_OOS.mat'),   'bootstrapqr');

switch cfg.model_selection
    case 2
        save(fullfile(sktDir,'skewtparam_inflation_OOS.mat'), ...
            'lc_skt','sc_skt','sh_skt','df_skt');
        save(fullfile(sktDir,'moments_inflation_OOS.mat'), ...
            'fcstmean','fcststdev','fcstskew');
    case 3
        save(fullfile(outDir,'semi_param','semi_param_distr_inflation.mat'), ...
            'semi_param_distr','-v7.3');
        save(fullfile(outDir,'semi_param','moments_inflation.mat'), ...
            'fcstmean','fcststdev','fcstskew');
        save(fullfile(outDir,'crps','crps_results_inflation.mat'), 'crps_results');
    case 4
        save(fullfile(outDir,'two_piece_normal','tpn_params_inflation.mat'), 'param_tpn');
        save(fullfile(outDir,'two_piece_normal','moments_inflation.mat'), ...
            'fcstmean','fcststdev','fcstskew');
end

% Export moments to Excel  (all model types)
xlsFile   = fullfile(outDir, sprintf('%s_best_spec_inflation.xlsx', ...
    cfg.modellist{cfg.model_selection}));
colNames  = cellstr(strcat('h_', string(0:cfg.horizons-1)));
for i = 1:3
    mList = {'fcstmean','fcststdev','fcstskew'};
    tbl   = [table(months_origin', 'VariableNames',{'Dates'}), ...
             array2table(eval(mList{i}), 'VariableNames', colNames)];
    writetable(tbl, xlsFile, 'Sheet', mList{i});
end

%% ════════════════════════════════════════════════════════════════════════
%%  7.  ROLLING FAN CHARTS  (one per horizon; x-axis = target date)
%% ════════════════════════════════════════════════════════════════════════

actualDT = datetime(dateNumeric_full, 'ConvertFrom','datenum');
qt_idx   = arrayfun(@(q) find(abs(cfg.quantiles-q)<1e-8,1), ...
                    [0.05 0.10 0.25 0.50 0.75 0.90 0.95]);

for h_plot = cfg.hPlot

    fcst_dates = months_origin + calmonths(h_plot - 1);
    Q = squeeze(pred_q(:, qt_idx, h_plot));   % nOrigins × 7

    [tf, loc] = ismember(datenum(fcst_dates), dateNumeric_full);
    act = NaN(size(fcst_dates));
    act(tf) = actual_var(loc(tf));

    fig = figure('Units','normalized','Position',[0.2 0.15 0.6 0.7],'Color','w');
    ax  = gca; hold(ax,'on');
    plotFanBands(ax, fcst_dates, Q);
    plot(ax, fcst_dates, Q(:,4), '--','Color',[0 0 0.7],'LineWidth',1.2,'DisplayName','50^{th}');
    plot(ax, fcst_dates, act,    'k', 'LineWidth',1.25, 'DisplayName','Outturn');
    styleAxis(ax, fcst_dates, 2);
    legend(ax,'show','Location','northoutside','Orientation','horizontal'); legend boxoff;
    set(ax,'FontSize',14);
    title(ax, sprintf('Inflation (CPI) — OOS fan  |  h = %d months ahead', h_plot-1));
    saveFig(fig, fanDir, sprintf('inflation_fanchart_OOS_h%d.png', h_plot));
end

%% ════════════════════════════════════════════════════════════════════════
%%  8.  LAST-ORIGIN FORWARD FAN  (1–36 months ahead)
%% ════════════════════════════════════════════════════════════════════════

t_last     = nOrigins;
lastOrigDT = months_origin(t_last);
n_fwd      = min(36, cfg.horizons - 1);
fcst_fwd   = lastOrigDT + calmonths(1:n_fwd);
Q_fwd      = squeeze(pred_q(t_last, qt_idx, 2:n_fwd+1))';   % n_fwd × 7

fig = figure('Units','normalized','Position',[0.18 0.12 0.64 0.72],'Color','w');
ax  = gca; hold(ax,'on');
plot(ax, actualDT, actual_var, 'k','LineWidth',1.25,'DisplayName','Outturn (CPI)');
xline(ax, lastOrigDT, 'k--','LineWidth',0.9,'DisplayName','Last origin');
plotFanBands(ax, fcst_fwd, Q_fwd);
plot(ax, fcst_fwd, Q_fwd(:,4), '--','Color',[0 0 0.7],'LineWidth',1.2,'DisplayName','Median (50^{th})');
yline(ax, 0,'k-','LineWidth',0.75,'HandleVisibility','off');
xlim(ax, [actualDT(max(end-120,1)), fcst_fwd(end)+calmonths(1)]);
styleAxis(ax, fcst_fwd, 2);
legend(ax,'show','Location','northoutside','Orientation','horizontal'); legend boxoff;
set(ax,'FontSize',13);
title(ax, sprintf('Inflation (CPI) — last-origin fan  |  1–%d months ahead', n_fwd));
saveFig(fig, fanDir, sprintf('inflation_fanchart_OOS_lastOrigin_%dmo.png', n_fwd));

%% ════════════════════════════════════════════════════════════════════════
%%  9.  HISTORICAL DECOMPOSITION
%% ════════════════════════════════════════════════════════════════════════

contrib    = NaN(nOrigins, Qn, numel(cfg.hPlot), V);
pred_check = NaN(nOrigins, Qn, numel(cfg.hPlot));

for t = 1:nOrigins
    for ih = 1:numel(cfg.hPlot)
        B = coeffqr(:,:, cfg.hPlot(ih), t);
        pred_check(t,:,ih) = Xfull(t,:) * B;
        for v = 1:V
            contrib(t,:,ih,v) = Xfull(t,v) * B(v,:);
        end
    end
end

q_idx_d = arrayfun(@(q) find(abs(cfg.quantiles-q)<1e-8,1), cfg.qDecomp);
cmap    = lines(V);

for qi = 1:numel(cfg.qDecomp)
    for ih = 1:numel(cfg.hPlot)
        h_plot = cfg.hPlot(ih);
        C      = squeeze(contrib(:, q_idx_d(qi), ih, :));
        qline  = squeeze(pred_check(:, q_idx_d(qi), ih));

        fig = figure('Units','normalized','Position',[0.2 0.15 0.6 0.7],'Color','w');
        ax  = gca; hold(ax,'on');
        b   = bar(ax, months_origin, C, 'stacked','BarWidth',0.7);
        for j = 1:V, b(j).FaceColor = cmap(j,:); end
        hl  = plot(ax, months_origin, qline, 'k-','LineWidth',1.3, ...
                   'DisplayName',sprintf('q=%.0f^{th}', cfg.qDecomp(qi)*100));
        yline(ax,0,'k-','LineWidth',0.75,'HandleVisibility','off');
        grid(ax,'on');
        ax.XTick = months_origin(1):calyears(2):months_origin(end);
        ax.XAxis.TickLabelFormat = 'yyyy';
        set(ax,'FontSize',13);
        legend(ax, [b(:)', hl], ...
               [vLabels, {sprintf('q=%.0f^{th}', cfg.qDecomp(qi)*100)}], ...
               'Location','southoutside','Orientation','horizontal','Box','off');
        title(ax, sprintf('Historical decomposition — inflation  |  h=%d, q=%.0f^{th}', ...
              h_plot-1, cfg.qDecomp(qi)*100));
        hold(ax,'off');
        saveFig(fig, figDir, sprintf('inflation_decomp_OOS_h%d_q%d.png', ...
                h_plot, round(cfg.qDecomp(qi)*100)));
    end
end

%% ════════════════════════════════════════════════════════════════════════
%%  10.  LÓPEZ-SALIDO CHART  (last vintage, bootstrap std-dev bands)
%%
%%   Layout:  5 predictors × 3 horizons  (1 month, 1 year, 2 years)
%%   Y-axis limits are set for the default specification; adjust if needed.
%% ════════════════════════════════════════════════════════════════════════

B_last  = squeeze(coeffqr(:,:,:,end));         % V × Qn × H
std_bst = squeeze(std(bootstrapqr, 0, 4));     % V × Qn × H

q_plot  = [0.05 0.10 0.25 0.50 0.75 0.90 0.95];
qt_ls   = arrayfun(@(q) find(abs(cfg.quantiles-q)<1e-8,1), q_plot);
qlabels = arrayfun(@(q) sprintf('%.0f^{th}',q*100), q_plot,'UniformOutput',false);
colors  = lines(numel(qt_ls));

% Predictor display order and labels  (indices into cfg.var.*)
var_order  = [1, 2, 3, 5, 4];   % persistence / expectations / slack / FCI / supply
var_labels = { sprintf('INFLATION\nPERSISTENCE'),  sprintf('INFLATION\nEXPECTATIONS'), ...
               sprintf('ECONOMIC\nTIGHTNESS'),     sprintf('FINANCIAL\nCONDITIONS'), ...
               sprintf('EXTERNAL\nCONDITIONS') };

% Y-axis limits (rows = var_order index; set to [] for auto-scale)
var_ylims  = { [-2,  1  ]; ...   % persistence
               [-1,  4  ]; ...   % expectations
               [-15, 15 ]; ...   % slack
               [-1.5e-2, 1.5e-2]; ... % FCI
               [-0.20, 0.40] };  % supply
var_yexp   = [0; 0; 0; -2; -2]; % YAxis.Exponent override (0 = no override)

horLabels  = {};
for iH = 1:numel(cfg.hor_ls)
    h = cfg.hor_ls(iH);
    if h == 2,     horLabels{end+1} = '1-MONTH';
    elseif h <= 4, horLabels{end+1} = sprintf('%d-MONTH', h-1);
    else,          horLabels{end+1} = sprintf('%g-YEAR', (h-1)/12);
    end
end

fig = figure('Units','centimeters','Position',[1 1 22 30],'Color','w');
tiledlayout(numel(var_order), numel(cfg.hor_ls), ...
    'TileSpacing','Compact','Padding','Compact');

for k = 1:numel(var_order)
    v_idx = 1 + var_order(k);    % +1 to skip constant row
    for iH = 1:numel(cfg.hor_ls)
        hor = cfg.hor_ls(iH);
        ax  = nexttile; hold(ax,'on'); yline(ax,0,'k');
        for dd = 1:numel(qt_ls)
            eb = errorbar(ax, dd, B_last(v_idx, qt_ls(dd), hor), ...
                          std_bst(v_idx, qt_ls(dd), hor), ...
                          's','LineWidth',1.5,'CapSize',0);
            eb.MarkerFaceColor = colors(dd,:);
            eb.MarkerEdgeColor = 'none';
        end
        ax.XLim       = [0.5, numel(qt_ls)+0.5];
        ax.XTick      = 1:numel(qt_ls);
        ax.XTickLabel = qlabels;
        if ~isempty(var_ylims{k}), ylim(ax, var_ylims{k}); end
        if var_yexp(k) ~= 0, ax.YAxis.Exponent = var_yexp(k); end
        ax.YAxis.TickLabelFormat = '%.1f';
        if iH == 1, ylabel(ax, var_labels{k}, 'Interpreter','none'); end
        if k  == 1
            title(ax, horLabels{iH}, 'Units','normalized','Position',[0.5,1.10,0]);
        end
        ax.FontSize = 10;  ax.Title.FontSize = 12;
        box(ax,'on');  hold(ax,'off');
    end
end
set(fig,'PaperUnits','centimeters','PaperPosition',[0 0 22 30]);
saveFig(fig, figDir, 'econ_interpr_inflation_OOS_1m_1y_2y.png', 300);

%% ════════════════════════════════════════════════════════════════════════
%%  11.  EXPORT LAST-ORIGIN SKEW-T PARAMETERS  (model == 2 only)
%% ════════════════════════════════════════════════════════════════════════

if cfg.model_selection == 2
    h_fwd_idx     = 2 : cfg.horizons;
    lc36 = lc_skt(t_last, h_fwd_idx)';   sc36 = sc_skt(t_last, h_fwd_idx)';
    sh36 = sh_skt(t_last, h_fwd_idx)';   df36 = df_skt(t_last, h_fwd_idx)';
    fcst_dates_36 = (lastOrigDT + calmonths(1:numel(h_fwd_idx)))';

    save(fullfile(sktDir,'last_origin_skewt_params_36mo.mat'), ...
        'lc36','sc36','sh36','df36','fcst_dates_36','lastOrigDT');

    % Quarterly cadence (every 3rd month → 12 quarters over 3 years)
    qidx = 3:3:numel(h_fwd_idx);
    save(fullfile(sktDir,'quarterly_skewt_params_forSharing.mat'), ...
        'lc36','sc36','sh36','df36','fcst_dates_36','lastOrigDT','-v7.3');
end

fprintf('\n── DONE: Inflation OOS pipeline complete. ──\n');

%% ════════════════════════════════════════════════════════════════════════
%%  LOCAL FUNCTIONS
%% ════════════════════════════════════════════════════════════════════════

function mkdirs(paths)
%MKDIRS  Create directories if they do not already exist.
    for i = 1:numel(paths)
        if ~exist(paths{i},'dir'), mkdir(paths{i}); end
    end
end

function plotFanBands(ax, dates, Q)
%PLOTFANBANDS  Draw shaded quantile fan bands.
%   dates : 1×T datetime;  Q : T×7 [Q05 Q10 Q25 Q50 Q75 Q90 Q95]
    d = dates(:)';
    bands = { [1,7], [0.85 0.9 1],  '5^{th}–95^{th}'; ...
              [2,6], [0.65 0.8 1],  '10^{th}–90^{th}'; ...
              [3,5], [0.4  0.6 1],  '25^{th}–75^{th}' };
    for b = 1:3
        lo = Q(:, bands{b,1}(1))';  hi = Q(:, bands{b,1}(2))';
        fill(ax, [d, fliplr(d)], [lo, fliplr(hi)], bands{b,2}, ...
             'EdgeColor','none','FaceAlpha',1,'DisplayName',bands{b,3});
    end
end

function styleAxis(ax, dates, step_yrs)
%STYLEAXIS  Zero line, grid, annual tick labels.
    yline(ax, 0,'k-','LineWidth',0.75,'HandleVisibility','off');
    grid(ax,'on');
    ticks = datetime(year(dates(1)),1,1):calyears(step_yrs):dates(end);
    ax.XTick = ticks;
    ax.XAxis.TickLabelFormat = 'yyyy';
end

function saveFig(fig, folder, filename, dpi)
%SAVEFIG  Save figure as PNG and close it.  Default dpi = 200.
    if nargin < 4, dpi = 200; end
    print(fig, fullfile(folder, filename), '-dpng', sprintf('-r%d',dpi));
    close(fig);
end
