%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  GDP AT-RISK  —  Recursive OOS Quantile LP + Model/Spec Selection
%%  (Quarterly, Multi-Specification, Covid dummy)
%
%  Authors : Aikman, Bidder, Lloyd, Mantoan, Maso, Mori, Tong
%  Updated : April 2026
%
%  PIPELINE
%    1. Generate all variable-category combinations  (combo_specifications)
%    2. For each spec: load data, estimate QR (with Covid dummy), build pred. quantiles
%    3. Select best spec via Weighted Interval Score (WIS)
%    4. For selected spec only: fit distribution (skew-t / semi-param / TPN)
%    5. Save .mat outputs
%    6. Figures: rolling fan / forward fans / HD / López-Salido (1Q,1Y,2Y)
%    7. Export last-origin skew-t parameters for sharing
%
%  MODEL SELECTION   (cfg.model_selection)
%    2 = QR + Skew-t (Azzalini-Capitanio)      [fitted at last origin only]
%    3 = QR + Semi-parametric (Mitchell-Poon-Zhu)
%    4 = QR + Two-piece Normal
%
%  SPEC SELECTION    (cfg.use_best_spec)
%    1 = auto: pick spec with lowest average WIS across all horizons
%    0 = manual: use cfg.specplot
%
%  Add alternatives to each cfg.var.* cell to search over more specs, e.g.:
%    cfg.var.current_act = {'pmi_out_long', 'mgdp_yoy'};
%
%  COVID DUMMY
%    Included only once the estimation window reaches cfg.covidDate.
%    The QR function drops the dummy column when collinear with the constant.
%
%  NOTE: Set cfg.bst.nboot = 5000 for production runs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% ── PATHS ────────────────────────────────────────────────────────────────
scriptDir = fileparts(mfilename('fullpath'));

dataFile = fullfile(scriptDir, 'GaRDataRaw_quarterly_BIS_march.xlsx');
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
cfg.startT   = datenum(1980, 6, 30);
cfg.endT     = datenum(2026, 3, 31);
cfg.startEst = datenum(2004, 3, 31);    % first forecast origin

% Dependent variable
cfg.var.dep = {'g4rgdp'};

% Variable categories — list alternatives in each cell to search over specs.
% The ndgrid below generates every combination automatically.
% e.g. two options in current_act × two in fci = 4 specs total.

cfg.var.current_act = {'mgdp_yoy'};               % current economic activity
%                      'pmi_out_long'              % (alternative: PMI output)
%                      'pmi_out_fut_long'          % (alternative: PMI output future)

cfg.var.leverage    = {'global_credit'};           % leverage / credit-to-GDP growth
%                      'delta_3y_credit_to_gdp_all'% (alternative: UK 3y credit growth)

cfg.var.fci         = {'ciss_uk'};                 % financial conditions
%                      'market_vol_uk'             % (alternative: market volatility)
%                      'yield_curve_slope'         % (alternative: yield curve slope)

cfg.var.macro_cond  = {'g4_import_deflator_fuel'}; % nominal / external indicator
%                      'g4infl'                    % (alternative: G4 inflation)
%                      'inflation_expectations'    % (alternative: infl expectations)

cfg.var.labour      = {'labour'};                  % labour-market indicator

% Distribution model
%   2 = Skew-t  |  3 = Semi-parametric  |  4 = Two-piece Normal
cfg.model_selection = 2;
cfg.modellist       = {'ols', 'skewt', 'semi-param', 'two-piece-normal'};

% Spec selection
cfg.use_best_spec = 1;       % 1 = auto (lowest avg WIS);  0 = manual
cfg.specplot      = 1;       % used when use_best_spec = 0
cfg.horizons_wis  = [1 5 9]; % horizons for per-horizon WIS table

% Covid dummy
cfg.covid     = 1;
cfg.covidDate = datenum(2020, 6, 30);   % first quarter the dummy is active

% Dependent-variable form
cfg.h_step_gdp = 0;  % 0 = year-on-year growth  |  1 = h-step annualised growth

% Model
cfg.horizons  = 13;            % h=1 current; h=2..13 → 1..12 quarters ahead
cfg.quantiles = 0.05:0.05:0.95;

% Bootstrap  *** SET nboot = 5000 FOR PRODUCTION RUNS ***
cfg.bst.nboot     = 10;
cfg.bst.blocksize = 8;
cfg.bst.ci        = 68;

% Plotting
cfg.hPlot   = [2, 5, 9];                         % horizons for rolling fan + HD
cfg.qDecomp = [0.05, 0.25, 0.50, 0.75, 0.90];   % quantiles for HD charts
cfg.hor_ls  = [2, 5, 9];                          % López-Salido (1Q, 1Y, 2Y)

%% ── GLOBAL PLOT DEFAULTS ─────────────────────────────────────────────────
set(0, 'defaultAxesFontName',      'Times');
set(0, 'defaultAxesLineStyleOrder', '-|--|:');
set(0, 'defaultLineLineWidth',       1);
rng('default');

%% ════════════════════════════════════════════════════════════════════════
%%  1.  GENERATE SPEC COMBINATIONS
%% ════════════════════════════════════════════════════════════════════════

[i0,i1,i2,i3,i4,i5] = ndgrid( ...
    1:numel(cfg.var.dep),         1:numel(cfg.var.current_act), ...
    1:numel(cfg.var.leverage),    1:numel(cfg.var.fci),         ...
    1:numel(cfg.var.macro_cond),  1:numel(cfg.var.labour));

vars = { cfg.var.dep(i0(:)),         cfg.var.current_act(i1(:)), ...
         cfg.var.leverage(i2(:)),    cfg.var.fci(i3(:)),         ...
         cfg.var.macro_cond(i4(:)),  cfg.var.labour(i5(:)) };
for k = 1:numel(vars)
    if size(vars{k}, 2) > 1, vars{k} = vars{k}.'; end
end

combo_specifications = [vars{:}];           % nSpec × (1 dep + 5 predictors)
nSpec = size(combo_specifications, 1);
nPred = size(combo_specifications, 2) - 1;  % predictors (excl. dep var)

% V_max: const + predictors + 1 Covid dummy row  (Covid row stays NaN pre-Covid)
V = nPred + 1 + 1;

spec_names = arrayfun(@(s) strjoin(combo_specifications(s,:), ' | '), ...
             (1:nSpec)', 'uni', false);
fprintf('Total specifications: %d\n', nSpec);

%% ════════════════════════════════════════════════════════════════════════
%%  2.  PRE-ALLOCATE ARRAYS
%% ════════════════════════════════════════════════════════════════════════

Qn   = numel(cfg.quantiles);
dv1  = datevec(cfg.startEst); dv2 = datevec(cfg.endT);
nOrigins = (dv2(1)-dv1(1))*4 + ceil(dv2(2)/3) - ceil(dv1(2)/3) + 1;

pred_q      = NaN(nOrigins, Qn, cfg.horizons, nSpec);
coeffqr     = NaN(V, Qn, cfg.horizons, nOrigins, nSpec);
bootstrapqr = NaN(V, Qn, cfg.horizons, cfg.bst.nboot, nSpec);
avg_wis     = NaN(nSpec, 1);

% Shared across specs (set once during first iteration)
idx_est          = [];
last_origin      = [];
idx_covid        = [];
dateNumeric_full = [];
actualvar        = [];

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

    dataloading_quarterly;

    %% ── First-iteration setup (spec-independent quantities) ─────────────
    if spec == 1
        idx_est          = find(dateNumeric == cfg.startEst, 1);
        last_origin      = numel(dateNumeric);
        idx_covid        = find(dateNumeric == cfg.covidDate, 1);
        dateNumeric_full = dateNumeric;
        assert(~isempty(idx_est), 'cfg.startEst not found in data.');

        % actualvar(t,h) = dep-var value at h-step-ahead target of origin t
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

    % Covid exogenous dummy (one-period lag to match BoE information set)
    exovar = [0; countryData.UK.covid(1:end-1)];

    %% ── Recursive QR loop ───────────────────────────────────────────────
    for endtime = idx_est : last_origin

        t = endtime - idx_est + 1;
        X = explvar(1:endtime, :);
        y = depvar( 1:endtime, :);

        % Build LP targets
        Y_LP = NaN(size(X,1), 1, cfg.horizons);
        if cfg.h_step_gdp == 0   % year-on-year growth
            for h = 1:cfg.horizons
                if size(X,1) > h
                    Y_LP(1:size(X,1)-h, 1, h) = y(1+h:end);
                end
            end
        else                     % h-step annualised growth
            for h = 1:cfg.horizons
                n_row = size(y,1) - h;
                if n_row > 0
                    Y_LP(1:n_row, 1, h) = 100 * (log(y(1+h:end)) - log(y(1:end-h))) * (4/h);
                end
            end
        end

        % Covid dummy logic
        if isempty(idx_covid) || endtime < idx_covid
            covid_now = 0;   exo_now = 0;
        else
            covid_now = 1;   exo_now = exovar(1:endtime, :);
        end

        Xmat        = [ones(size(X,1),1), X];
        doBootstrap = (endtime == last_origin);

        [bQR, bQRbst] = qfe_qr_local_projection_SL_final( ...
            Y_LP, Xmat, cfg.quantiles, (1:cfg.horizons)', ...
            doBootstrap, cfg.bst, covid_now, exo_now);

        nCoeff = size(bQR, 1);
        coeffqr(1:nCoeff,:,:,t,spec) = bQR;
        if doBootstrap, bootstrapqr(1:nCoeff,:,:,:,spec) = bQRbst; end

        % Predicted quantiles at this origin
        if covid_now == 0
            x_now = [1, X(endtime,:)];
        else
            x_now = [1, X(endtime,:), exovar(endtime)];
        end
        x_padded = zeros(1, V);
        x_padded(1:nCoeff) = x_now(1:nCoeff);

        for h = 1:cfg.horizons
            pred_q(t,:,h,spec) = x_padded(1:nCoeff) * bQR(:,:,h);
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
        pq_h = squeeze(pred_q(:,:,h,s));
        av_h = actualvar(:,h);
        mask = ~isnan(av_h);
        if any(mask)
            u_h = repmat(av_h(mask), [1,Qn]) - pq_h(mask,:);
            ql_ = repmat(cfg.quantiles, [sum(mask),1]);
            wis_by_hor(s,h) = mean(max(ql_.*u_h, (ql_-1).*u_h), 'all');
        end
    end
end
[~, idx_min_hor] = min(wis_by_hor(:, cfg.horizons_wis), [], 1);

% Save WIS tables to Excel
xlsWIS = fullfile(wisDir, sprintf('wis_%s_gdp.xlsx', cfg.modellist{cfg.model_selection}));
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
dataloading_quarterly;

vlag    = strcat('l1', varnames);
depvar  = countryData.UK.(vlag{1});
explvar = table2array(countryData.UK(:, vlag(2:end)));
exovar  = [0; countryData.UK.covid(1:end-1)];
actual_var = countryData.UK.(cfg.var.dep{1});

StartEstDT     = datetime(cfg.startEst,              'ConvertFrom','datenum');
lastOriDT      = datetime(dateNumeric_full(last_origin), 'ConvertFrom','datenum');
quarters_origin = StartEstDT : calquarters(1) : lastOriDT;  % 1 × nOrigins

Xorig  = explvar(idx_est:last_origin, :);
exorig = exovar( idx_est:last_origin, :);
Xfull  = [ones(nOrigins,1), Xorig, exorig];   % nOrigins × V
vLabels = [{'Constant'}, ...
    cellfun(@(s) strrep(s,'_',' '), combo_specifications(spec_to_use,2:end), 'uni',0), ...
    {'Covid dummy'}];

save(fullfile(outDir,'explanatoryvar_gdp_OOS.mat'), 'explvar', 'exovar');

%% ════════════════════════════════════════════════════════════════════════
%%  5.  DISTRIBUTION FITTING  (for spec_to_use; skew-t at last origin only)
%% ════════════════════════════════════════════════════════════════════════

t_last = nOrigins;

switch cfg.model_selection

    case 2  %────────────────── Skew-t (last origin only, for speed) ──────
        PQ = reshape(pred_q(t_last,:,:), 1, 1, Qn, cfg.horizons, 1);
        [aa, ~, ~, ~, bb] = fit_skewt_to_quantiles_all( ...
            PQ, cfg.quantiles, [0.05 0.25 0.5 0.75 0.95], 1, 1:cfg.horizons);
        param_skt = squeeze(aa);               % 4 × H  [lc; sc; sh; df]
        fcstmean  = squeeze(bb(:,1,:,:))';     % H × 1
        fcststdev = sqrt(squeeze(bb(:,2,:,:))');
        fcstskew  = squeeze(bb(:,3,:,:))';

    case 3  %────────────────── Semi-parametric ───────────────────────────
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

save(fullfile(outDir,'actual_gdp_yoy_OOS.mat'), ...
    'actual_var','dateNumeric_full','idx_est','last_origin','startT','endT');
save(fullfile(outDir,'qreg_results_gdp_OOS.mat'),        'coeffqr');
save(fullfile(outDir,'predicted_quantiles_gdp_OOS.mat'), 'pred_q');
save(fullfile(outDir,'bootstrap_results_gdp_OOS.mat'),   'bootstrapqr');

switch cfg.model_selection
    case 2
        save(fullfile(sktDir,'sktparam_gdp_OOS.mat'), 'param_skt');
    case 3
        save(fullfile(outDir,'semi_param','moments_gdp.mat'), ...
            'fcstmean','fcststdev','fcstskew');
        save(fullfile(outDir,'crps','crps_results_gdp.mat'), 'crps_results');
    case 4
        save(fullfile(outDir,'two_piece_normal','tpn_params_gdp.mat'), 'param_tpn');
end

% Export moments to Excel
xlsFile  = fullfile(outDir, sprintf('%s_best_spec_gdp.xlsx', ...
    cfg.modellist{cfg.model_selection}));
colNames = cellstr(strcat('h_', string(0:cfg.horizons-1)));
momentVars = {'fcstmean','fcststdev','fcstskew'};
for i = 1:3
    tbl = [table(quarters_origin', 'VariableNames',{'Dates'}), ...
           array2table(eval(momentVars{i}), 'VariableNames', colNames)];
    writetable(tbl, xlsFile, 'Sheet', momentVars{i});
end

%% ════════════════════════════════════════════════════════════════════════
%%  7.  ROLLING FAN CHARTS  (one per horizon; x-axis = target date)
%% ════════════════════════════════════════════════════════════════════════

actualDT = datetime(dateNumeric_full, 'ConvertFrom','datenum');
qt_idx   = arrayfun(@(q) find(abs(cfg.quantiles-q)<1e-8,1), ...
                    [0.05 0.10 0.25 0.50 0.75 0.90 0.95]);

for h_plot = cfg.hPlot

    fcst_dates = quarters_origin + calquarters(h_plot - 1);
    Q = squeeze(pred_q(:, qt_idx, h_plot));

    [tf, loc] = ismember(datenum(fcst_dates), dateNumeric_full);
    act = NaN(size(fcst_dates));
    act(tf) = actual_var(loc(tf));

    fig = figure('Units','normalized','Position',[0.2 0.15 0.6 0.7],'Color','w');
    ax  = gca; hold(ax,'on');
    plotFanBands(ax, fcst_dates, Q);
    plot(ax, fcst_dates, Q(:,4), '--','Color',[0 0 0.7],'LineWidth',1.2,'DisplayName','50^{th}');
    plot(ax, fcst_dates, act,    'k', 'LineWidth',1.25,'DisplayName','Outturn');
    styleAxis(ax, fcst_dates, 2);
    legend(ax,'show','Location','northoutside','Orientation','horizontal'); legend boxoff;
    set(ax,'FontSize',14);
    title(ax, sprintf('GDP (YoY) — OOS fan  |  h = %d quarters ahead', h_plot-1));
    saveFig(fig, fanDir, sprintf('gdp_fanchart_OOS_h%d.png', h_plot));
end

%% ════════════════════════════════════════════════════════════════════════
%%  8.  LAST-ORIGIN FORWARD FANS  (1-year and 3-year)
%% ════════════════════════════════════════════════════════════════════════

lastOrigDT = quarters_origin(t_last);

for n_fwd = [4, 12]

    fcst_fwd = lastOrigDT + calquarters(1:n_fwd);
    Q_fwd    = squeeze(pred_q(t_last, qt_idx, 2:n_fwd+1))';   % n_fwd × 7

    fig = figure('Units','normalized','Position',[0.18 0.12 0.64 0.72],'Color','w');
    ax  = gca; hold(ax,'on');
    plot(ax, actualDT, actual_var, 'k','LineWidth',1.25,'DisplayName','Outturn (GDP YoY)');
    xline(ax, lastOrigDT, 'k--','LineWidth',0.9,'DisplayName','Last origin');
    plotFanBands(ax, fcst_fwd, Q_fwd);
    plot(ax, fcst_fwd, Q_fwd(:,4), '--','Color',[0 0 0.7],'LineWidth',1.2,'DisplayName','Median (50^{th})');
    yline(ax, 0,'k-','LineWidth',0.75,'HandleVisibility','off');
    xlim(ax, [datetime(2015,1,1), fcst_fwd(end)+calquarters(1)]);
    styleAxis(ax, fcst_fwd, 2);
    legend(ax,'show','Location','northoutside','Orientation','horizontal'); legend boxoff;
    set(ax,'FontSize',13);
    title(ax, sprintf('GDP (YoY) — last-origin fan  |  1–%d quarters ahead', n_fwd));
    saveFig(fig, fanDir, sprintf('gdp_fanchart_OOS_lastOrigin_%dQ.png', n_fwd));
end

%% ════════════════════════════════════════════════════════════════════════
%%  9.  HISTORICAL DECOMPOSITION
%%
%%  Pre-Covid origins: Covid coefficient is NaN → replace with 0 so that
%%  the dummy's contribution is zero (not NaN) for those periods.
%% ════════════════════════════════════════════════════════════════════════

contrib    = NaN(nOrigins, Qn, numel(cfg.hPlot), V);
pred_check = NaN(nOrigins, Qn, numel(cfg.hPlot));

for t = 1:nOrigins
    for ih = 1:numel(cfg.hPlot)
        B = coeffqr(:,:, cfg.hPlot(ih), t);
        B(isnan(B)) = 0;    % pre-Covid origins: treat missing Covid coeff as zero
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
        b   = bar(ax, quarters_origin, C, 'stacked','BarWidth',0.7);
        for j = 1:V, b(j).FaceColor = cmap(j,:); end
        hl  = plot(ax, quarters_origin, qline, 'k-','LineWidth',1.3, ...
                   'DisplayName',sprintf('q=%.0f^{th}', cfg.qDecomp(qi)*100));
        yline(ax,0,'k-','LineWidth',0.75,'HandleVisibility','off');
        grid(ax,'on');
        ax.XTick = quarters_origin(1):calyears(2):quarters_origin(end);
        ax.XAxis.TickLabelFormat = 'yyyy';
        set(ax,'FontSize',13);
        legend(ax, [b(:)', hl], ...
               [vLabels, {sprintf('q=%.0f^{th}', cfg.qDecomp(qi)*100)}], ...
               'Location','southoutside','Orientation','horizontal','Box','off');
        title(ax, sprintf('Historical decomposition — GDP  |  h=%d, q=%.0f^{th}', ...
              h_plot-1, cfg.qDecomp(qi)*100));
        hold(ax,'off');
        saveFig(fig, figDir, sprintf('gdp_decomp_OOS_h%d_q%d.png', ...
                h_plot, round(cfg.qDecomp(qi)*100)));
    end
end

%% ════════════════════════════════════════════════════════════════════════
%%  10.  LÓPEZ-SALIDO CHART  (last vintage, bootstrap std-dev bands)
%%
%%   Layout:  5 predictors + Covid dummy  ×  3 horizons  (1Q, 1Y, 2Y)
%%   Y-axis limits are set for the default specification; adjust if needed.
%% ════════════════════════════════════════════════════════════════════════

B_last  = squeeze(coeffqr(:,:,:,end));         % V × Qn × H
std_bst = squeeze(std(bootstrapqr, 0, 4));     % V × Qn × H

q_plot  = [0.05 0.10 0.25 0.50 0.75 0.90 0.95];
qt_ls   = arrayfun(@(q) find(abs(cfg.quantiles-q)<1e-8,1), q_plot);
qlabels = arrayfun(@(q) sprintf('%.0f^{th}',q*100), q_plot,'UniformOutput',false);
colors  = lines(numel(qt_ls));

% Predictor display order, labels and fixed y-limits
% Row order: economic activity, leverage, FCI, macro/nominal, labour, Covid
var_order  = 1 : nPred + 1;          % all predictors + Covid dummy
var_labels = { sprintf('ECONOMIC\nACTIVITY'),   sprintf('LEVERAGE\nGROWTH'),   ...
               sprintf('FINANCIAL\nCONDITIONS'), sprintf('NOMINAL\nINDICATOR'), ...
               sprintf('LABOUR\nMARKET'),         sprintf('COVID\nDUMMY') };
var_ylims  = { [-100, 100]; [-0.6, 0.2]; [-25, 20]; [-0.1, 0.1]; []; [] };
%              econ act       leverage     FCI          macro/nom  labour  Covid
%              (set to [] for auto-scale)

horLabels = {};
for iH = 1:numel(cfg.hor_ls)
    h = cfg.hor_ls(iH);
    if h == 2, horLabels{end+1} = '1-QUARTER';
    else,      horLabels{end+1} = sprintf('%g-YEAR', (h-1)/4);
    end
end

fig = figure('Units','centimeters','Position',[1 1 22 28],'Color','w');
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
        if k <= numel(var_ylims) && ~isempty(var_ylims{k}), ylim(ax, var_ylims{k}); end
        ax.YAxis.TickLabelFormat = '%.1f';
        if iH == 1, ylabel(ax, var_labels{k}, 'Interpreter','none'); end
        if k  == 1
            title(ax, horLabels{iH}, 'Units','normalized','Position',[0.5,1.10,0]);
        end
        ax.FontSize = 10;  ax.Title.FontSize = 12;
        box(ax,'on');  hold(ax,'off');
    end
end
set(fig,'PaperUnits','centimeters','PaperPosition',[0 0 22 28]);
saveFig(fig, figDir, 'econ_interpr_gdp_OOS_1q_1y_2y.png', 300);

%% ════════════════════════════════════════════════════════════════════════
%%  11.  EXPORT LAST-ORIGIN SKEW-T PARAMETERS  (model == 2 only)
%% ════════════════════════════════════════════════════════════════════════

if cfg.model_selection == 2
    h_fwd_idx  = 2 : cfg.horizons;
    lc_qtr = param_skt(1, h_fwd_idx)';   sc_qtr = param_skt(2, h_fwd_idx)';
    sh_qtr = param_skt(3, h_fwd_idx)';   df_qtr = param_skt(4, h_fwd_idx)';
    fcst_dates_qtr = (lastOrigDT + calquarters(1:numel(h_fwd_idx)))';

    save(fullfile(sktDir,'last_origin_skewt_params_12Q.mat'), ...
        'lc_qtr','sc_qtr','sh_qtr','df_qtr', ...
        'fcstmean','fcststdev','fcstskew','fcst_dates_qtr','lastOrigDT');
    save(fullfile(sktDir,'quarterly_skewt_params_forSharing_gdp.mat'), ...
        'lc_qtr','sc_qtr','sh_qtr','df_qtr','fcst_dates_qtr','lastOrigDT','-v7.3');
end

fprintf('\n── DONE: GDP OOS pipeline complete. ──\n');

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
