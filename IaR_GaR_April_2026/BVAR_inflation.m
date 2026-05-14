%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  INFLATION AT-RISK  —  Recursive OOS BVAR  (Monthly, Minnesota Prior)
%%
%%  Authors : Aikman, Bidder, Lloyd, Mantoan, Maso, Mori, Tong
%%  Updated : April 2026
%%
%%  PIPELINE
%%    1. Load raw VAR variables for the baseline specification
%%    2. Estimate prior scale factors (sigma) from AR(1) residuals
%%    3. Recursive OOS: for each origin, estimate on expanding window,
%%       build Minnesota prior via dummy observations, compute Normal-IW
%%       posterior, draw from posterior, simulate h-step paths
%%    4. Collect predictive quantiles  →  pred_q  (same format as QR)
%%    5. Compute average WIS for comparison against QR specifications
%%    6. Save .mat output; produce rolling and forward fan charts
%%
%%  BVAR SETTINGS  (cfg.bvar.*)
%%    p       : VAR lag order
%%    lambda1 : overall tightness  (smaller = tighter; 0.1–0.3 typical)
%%    lambda3 : lag decay exponent (1 = harmonic; 2 = faster decay)
%%    lambda4 : diffuseness of constant prior (large = very diffuse)
%%    delta   : N×1 prior mean for own first lag (1 = RW; 0 = stationary)
%%    nDraws  : posterior Monte Carlo draws per origin
%%
%%  OUTPUT FORMAT
%%    pred_q  : nOrigins × Qn × horizons  — same as QR scripts
%%    Variable 1 in the VAR (cfg.var.dep) is the forecast target.
%%    h=1 = nowcast; h=2..37 = 1..36 months ahead.
%%
%%  REQUIRES: Statistics and Machine Learning Toolbox (wishrnd)
%%            dataloading_monthly.m  (intermediate_codes folder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%% ── PATHS ────────────────────────────────────────────────────────────────
scriptDir = fileparts(mfilename('fullpath'));

dataFile = fullfile(scriptDir, 'IaRDataRaw_monthly_M.xlsx');
outDir   = fullfile(scriptDir, 'Outputs', 'BVAR');
figDir   = fullfile(outDir, 'figures_inflation');

%mkdirs({outDir, figDir});

addpath(fullfile(scriptDir, 'intermediate_codes'));
addpath(fullfile(scriptDir, 'functions'));

%% ── CONFIGURATION ────────────────────────────────────────────────────────

% Sample
cfg.startT   = datenum(1990, 2, 1);
cfg.endT     = datenum(2026, 3, 1);
cfg.startEst = datenum(2004, 1, 1);    % first forecast origin

% VAR variables — variable 1 is the forecast target (dep var).
% Swap in commented alternatives to change the specification.
cfg.var.dep     = 'g4cpi';

cfg.var.persist = 'avgcpi';                       % inflation persistence
%                  'avg12Infl'

cfg.var.expect  = 'infl1_m';                      % inflation expectations
%                  'infl2_m'

cfg.var.slack   = 'delta1vu_tightness';           % economic slack
%                  'ugap_hp_filter_lambda_129600'
%                  'ugap_kalman_filter'

cfg.var.supply  = 'yoy_growth_import_deflator';   % supply / external
%                  'g4oil'
%                  'g4cpicore_global'

cfg.var.fci     = 'bond_spread';                  % financial conditions
%                  'market_vol_uk'

% Collect in VAR order (dep var first, then predictors)
cfg.varnames = {cfg.var.dep, cfg.var.persist, cfg.var.expect, ...
                cfg.var.slack, cfg.var.supply, cfg.var.fci};

% Minnesota prior hyperparameters
cfg.bvar.p       = 1;     % VAR lag order  (1 recommended for monthly data)
cfg.bvar.lambda1 = 0.2;   % overall tightness
cfg.bvar.lambda3 = 1;     % lag decay exponent (1 = harmonic decay)
cfg.bvar.lambda4 = 1e3;   % constant prior scale (large = diffuse)
cfg.bvar.delta   = ones(numel(cfg.varnames), 1);  % own-lag prior mean
%   Set delta(n) = 0 for variables already in growth-rate / stationary form.
%   Set delta(n) = 1 for levels / near-unit-root series.

% Posterior Monte Carlo
cfg.bvar.nDraws = 2000;   % draws per origin  *** SET >= 5000 FOR PRODUCTION ***

% Forecast settings
cfg.horizons  = 37;             % h=1 nowcast; h=2..37 → 1..36 months ahead
cfg.quantiles = 0.05:0.05:0.95;

% Plotting
cfg.hPlot = [2, 4, 13];   % horizons shown in rolling fan and forward fan

%% ── GLOBAL PLOT DEFAULTS ─────────────────────────────────────────────────
set(0, 'defaultAxesFontName',      'Times');
set(0, 'defaultAxesLineStyleOrder', '-|--|:');
set(0, 'defaultLineLineWidth',       1);
rng('default');

%% ════════════════════════════════════════════════════════════════════════
%%  1.  LOAD DATA  (no extra lag columns — VAR handles its own lags)
%% ════════════════════════════════════════════════════════════════════════

N = numel(cfg.varnames);

varnames     = cfg.varnames;
lagctryvar   = zeros(1, N);   % no extra lag columns in the table
lagglobalvar = [];
lags         = lagctryvar;
startT       = cfg.startT;   endT         = cfg.endT;
ctrynames    = {'UK'};        onlyuk       = 1;
fullFileName = dataFile;      covid        = 0;
dummyvarname = {'covid'};

dataloading_monthly;   % produces countryData and dateNumeric

% Full-sample data matrix  (T_full × N)
T_full = height(countryData.UK);
Y_all  = nan(T_full, N);
for n = 1:N
    Y_all(:, n) = countryData.UK.(cfg.varnames{n});
end

%% ── Estimation start and origin count ───────────────────────────────────
idx_est     = find(dateNumeric == cfg.startEst, 1);
last_origin = T_full;
assert(~isempty(idx_est), 'cfg.startEst not found in dateNumeric.');

nOrigins = last_origin - idx_est + 1;

%% ── Actual dep-var matrix for WIS  (nOrigins × horizons) ────────────────
act_filt = Y_all(idx_est:end, 1);
n_act    = numel(act_filt);
interm   = NaN(n_act, n_act);
for ii = 1:n_act
    interm(1:(n_act-ii+1), ii) = act_filt(ii:end);
end
actualvar = interm(1:cfg.horizons, :)';   % nOrigins × horizons

%% ════════════════════════════════════════════════════════════════════════
%%  2.  PRIOR SCALE FACTORS  (sigma from AR(1) residual SDs)
%% ════════════════════════════════════════════════════════════════════════
%  Estimated once on the full sample.  These scale the Minnesota prior
%  so that shrinkage is proportional to each variable's own variability.

sigma_ar = nan(N, 1);
for n = 1:N
    y_n = Y_all(:, n);
    y_n = y_n(~isnan(y_n));
    if numel(y_n) < 4
        sigma_ar(n) = std(y_n, 'omitnan');
        continue
    end
    Z_ar      = [ones(numel(y_n)-1, 1), y_n(1:end-1)];
    b_ar      = Z_ar \ y_n(2:end);
    sigma_ar(n) = std(y_n(2:end) - Z_ar * b_ar);
end

%% ════════════════════════════════════════════════════════════════════════
%%  3.  PRE-ALLOCATE OUTPUT
%% ════════════════════════════════════════════════════════════════════════

Qn     = numel(cfg.quantiles);
p      = cfg.bvar.p;
K      = 1 + N * p;           % regressors per equation: const + N*p lags

pred_q = NaN(nOrigins, Qn, cfg.horizons);

%% ════════════════════════════════════════════════════════════════════════
%%  4.  RECURSIVE OOS LOOP
%% ════════════════════════════════════════════════════════════════════════
%
%  At origin o  (absolute row = idx_est + o - 1):
%    Information set : Y_all(1 : origin_abs - 1, :)   [data through t-1]
%    Forecast target : column 1 of Y_all               [dep var = cfg.var.dep]
%    h=1             : nowcast (period t = origin_abs)
%    h=2..H          : 1..H-1 months ahead
%
%  Steps for each origin:
%    (a) Build Y_est, X_est from the estimation window
%    (b) Construct Minnesota dummy observations (Yd, Xd)
%    (c) OLS on augmented system → posterior mode B_hat and scale S_post
%    (d) Draw (Sigma, B) from Normal-IW posterior
%    (e) Simulate H-step forward path; record dep-var forecasts
%    (f) Compute quantiles across draws → pred_q(o, :, h)

fprintf('BVAR inflation: running %d origins...\n', nOrigins);
tic;

for o = 1:nOrigins

    origin_abs = idx_est + o - 1;     % current origin (absolute row)
    T_full_est = origin_abs - 1;      % data available: rows 1..T_full_est

    %% ── Build clean estimation sample (drop NaN rows) ────────────────────
    Y_raw = Y_all(1:T_full_est, :);
    ok    = all(~isnan(Y_raw), 2);
    Y_est = Y_raw(ok, :);
    T_est = size(Y_est, 1);

    if T_est <= p + K + 2, continue; end   % skip if too few obs

    %% ── Construct LHS and RHS matrices ──────────────────────────────────
    T_reg = T_est - p;
    Y_lhs = Y_est(p+1:end, :);        % T_reg × N
    X_rhs = zeros(T_reg, K);
    X_rhs(:, 1) = 1;                  % constant
    for lag = 1:p
        cols = 1 + (lag-1)*N + (1:N);
        X_rhs(:, cols) = Y_est(p+1-lag : T_est-lag, :);
    end

    %% ── Minnesota dummy observations ─────────────────────────────────────
    %
    %  Dummy 1  (N*p rows): prior on lag coefficients
    %    Row (l-1)*N+j: Yd = sigma(j)*delta(j) / (l^lambda3 * lambda1)  [col j]
    %                   Xd = sigma(j)           / (l^lambda3 * lambda1)  [lag-l block, col j]
    %
    %  Dummy 2  (N rows): prior on error covariance  (Yd = diag(sigma))
    %
    %  Dummy 3  (1 row):  diffuse prior on constant  (Xd(1,1) = 1/lambda4)

    lam1 = cfg.bvar.lambda1;
    lam3 = cfg.bvar.lambda3;
    lam4 = cfg.bvar.lambda4;
    delt = cfg.bvar.delta;
    sig  = sigma_ar;

    % Dummy 1
    Yd1 = zeros(N*p, N);
    Xd1 = zeros(N*p, K);
    for l = 1:p
        for j = 1:N
            row = (l-1)*N + j;
            scale = sig(j) / (l^lam3 * lam1);
            Yd1(row, j)                  = scale * delt(j);
            Xd1(row, 1 + (l-1)*N + j)   = scale;
        end
    end

    % Dummy 2
    Yd2 = diag(sig);
    Xd2 = zeros(N, K);

    % Dummy 3
    Yd3    = zeros(1, N);
    Xd3    = zeros(1, K);
    Xd3(1) = 1 / lam4;

    % Augmented system
    Ystar = [Yd1; Yd2; Yd3; Y_lhs];
    Xstar = [Xd1; Xd2; Xd3; X_rhs];

    %% ── Normal-IW posterior (OLS on augmented data) ──────────────────────
    %
    %  Posterior mean   : B_hat = (Xstar'Xstar)^{-1} Xstar'Ystar   [K × N]
    %  Posterior scale  : S_post = (Ystar - Xstar B_hat)'(...)      [N × N]
    %  Posterior df     : nu_post = rows(Ystar) - K

    XXi    = (Xstar' * Xstar) \ eye(K);         % K×K  (posterior Omega)
    B_hat  = XXi * (Xstar' * Ystar);            % K×N  posterior mean
    E_aug  = Ystar - Xstar * B_hat;
    S_post = E_aug' * E_aug;                    % N×N  posterior scale
    nu_post = size(Ystar, 1) - K;

    if nu_post <= N + 1, continue; end           % degenerate — skip

    %% ── Pre-compute factors for efficient draws ──────────────────────────
    [L_Omega, flag_Om] = chol(XXi, 'lower');    % K×K  chol of posterior covariance
    if flag_Om > 0, continue; end               % XXi not PD — skip origin
    S_inv = S_post \ eye(N);                    % for wishrnd

    %% ── Last p rows of Y_est as forecast initial conditions ──────────────
    Y_init = Y_est(end-p+1:end, :);             % p×N

    %% ── Monte Carlo predictive draws ─────────────────────────────────────
    %
    %  Numerically stable IW draw — avoids explicit inversion of W_draw:
    %    W_draw ~ Wishart(S_inv, nu_post)
    %    W_draw = L_W * L_W'  (lower Cholesky)
    %    U = inv(L_W')        (upper triangular)
    %    Then  U * U' = Sigma_d = inv(W_draw)  [no inv() needed]
    %    Draws from N(0, Sigma_d): U * randn(N,1)
    %    B draw: L_Omega * randn(K,N) * U'

    H      = cfg.horizons;
    nD     = cfg.bvar.nDraws;
    fc_dep = NaN(nD, H);

    for d = 1:nD

        % Draw W ~ Wishart(S_inv, nu_post); get its Cholesky
        W_draw        = wishrnd(S_inv, nu_post);
        [L_W, flag_W] = chol(W_draw, 'lower');
        if flag_W > 0, continue; end        % non-PD draw — discard

        % U (upper triangular) satisfies U*U' = Sigma_d = inv(W_draw)
        U = L_W' \ eye(N);                  % stable triangular solve

        % Draw B | Sigma ~ MN(B_hat, Sigma_d ⊗ XXi)
        B_d = B_hat + L_Omega * randn(K, N) * U';

        % Simulate H steps forward
        Y_path = [Y_init; NaN(H, N)];       % (p+H) × N
        for h = 1:H
            z_h    = zeros(1, K);
            z_h(1) = 1;
            for l = 1:p
                cols_l = 1 + (l-1)*N + (1:N);
                z_h(cols_l) = Y_path(p + h - l, :);
            end
            % Shock: U * z ~ N(0, Sigma_d)
            Y_path(p + h, :) = z_h * B_d + (U * randn(N, 1))';
        end

        fc_dep(d, :) = Y_path(p+1 : p+H, 1)';
    end

    %% ── Store predictive quantiles ───────────────────────────────────────
    for h = 1:H
        fc_h = fc_dep(:, h);
        fc_h = fc_h(~isnan(fc_h));
        if numel(fc_h) < 10, continue; end
        pred_q(o, :, h) = quantile(fc_h, cfg.quantiles);
    end

    if mod(o, 20) == 0
        fprintf('  Origin %3d / %d  (%.0f s)\n', o, nOrigins, toc);
    end

end

fprintf('Done. Total time: %.1f min\n', toc / 60);

%% ════════════════════════════════════════════════════════════════════════
%%  5.  WEIGHTED INTERVAL SCORE  (for comparison with QR specifications)
%% ════════════════════════════════════════════════════════════════════════
ql3      = reshape(cfg.quantiles, 1, Qn, 1);
av3      = repmat(reshape(actualvar, nOrigins, 1, cfg.horizons), [1, Qn, 1]);
u3       = av3 - pred_q;
check    = max(ql3 .* u3, (ql3 - 1) .* u3);
bvar_wis = mean(check(~isnan(av3(:))), 'omitnan');
fprintf('BVAR average WIS: %.4f\n', bvar_wis);

%% ════════════════════════════════════════════════════════════════════════
%%  6.  SAVE
%% ════════════════════════════════════════════════════════════════════════
dateNumeric_est = dateNumeric(idx_est : idx_est + nOrigins - 1);

save(fullfile(outDir, 'BVAR_inflation_pred_q.mat'), ...
    'pred_q', 'actualvar', 'dateNumeric_est', 'cfg', 'bvar_wis', 'sigma_ar');

fprintf('Saved: %s\n', fullfile(outDir, 'BVAR_inflation_pred_q.mat'));

%% ════════════════════════════════════════════════════════════════════════
%%  7.  FIGURES
%% ════════════════════════════════════════════════════════════════════════

qt_idx     = arrayfun(@(q) find(abs(cfg.quantiles-q)<1e-8,1), ...
                      [0.05 0.10 0.25 0.50 0.75 0.90 0.95]);

months_est = datetime(dateNumeric_est, 'ConvertFrom','datenum');
actualDT   = datetime(dateNumeric,     'ConvertFrom','datenum');
actual_var = Y_all(:, 1);   % full dep-var series

%% ── (a) Rolling fan charts  (x-axis = target date) ──────────────────────
for ih = 1:numel(cfg.hPlot)
    h = cfg.hPlot(ih);
    Q = squeeze(pred_q(:, qt_idx, h));   % nOrigins × 7

    tgt_dates = months_est + calmonths(h - 1);

    valid  = ~all(isnan(Q), 2);
    Q_v    = Q(valid, :);
    tgt_v  = tgt_dates(valid);

    oos_mask = actualDT >= tgt_dates(1);

    fig = figure('Units','normalized','Position',[0.2 0.15 0.6 0.7],'Color','w');
    ax  = gca; hold(ax,'on');
    plotFanBands(ax, tgt_v, Q_v);
    plot(ax, tgt_v, Q_v(:,4), '--', 'Color',[0 0 0.7], 'LineWidth',1.2, ...
         'DisplayName','50^{th}');
    plot(ax, actualDT(oos_mask), actual_var(oos_mask), 'k', 'LineWidth',1.25, ...
         'DisplayName','Outturn');
    styleAxis(ax, tgt_dates, 2);
    legend(ax,'show','Location','northoutside','Orientation','horizontal'); legend boxoff;
    set(ax,'FontSize',14);
    title(ax, sprintf('Inflation (CPI) — BVAR rolling fan  |  h = %d months ahead', h-1));
    saveFig(fig, figDir, sprintf('BVAR_infl_rolling_h%02d.png', h));
end

%% ── (b) Forward predictive density  (last available origin) ─────────────
lastOrigDT = months_est(end);
n_fwd      = min(36, cfg.horizons - 1);
fcst_fwd   = lastOrigDT + calmonths(1:n_fwd);
Q_fwd      = squeeze(pred_q(nOrigins, qt_idx, 2:n_fwd+1))';   % n_fwd × 7

fig = figure('Units','normalized','Position',[0.18 0.12 0.64 0.72],'Color','w');
ax  = gca; hold(ax,'on');
plot(ax, actualDT, actual_var, 'k','LineWidth',1.25,'DisplayName','Outturn (CPI)');
xline(ax, lastOrigDT, 'k--','LineWidth',0.9,'DisplayName','Last origin');
plotFanBands(ax, fcst_fwd, Q_fwd);
plot(ax, fcst_fwd, Q_fwd(:,4), '--','Color',[0 0 0.7],'LineWidth',1.2, ...
     'DisplayName','Median (50^{th})');
yline(ax, 0,'k-','LineWidth',0.75,'HandleVisibility','off');
xlim(ax, [actualDT(max(end-120,1)), fcst_fwd(end)+calmonths(1)]);
styleAxis(ax, fcst_fwd, 2);
legend(ax,'show','Location','northoutside','Orientation','horizontal'); legend boxoff;
set(ax,'FontSize',13);
title(ax, sprintf('Inflation (CPI) — BVAR forward fan  |  1–%d months ahead', n_fwd));
saveFig(fig, figDir, 'BVAR_infl_forward_fan.png');

%% ── (c) Median predictive heatmap ───────────────────────────────────────
fig = figure('Visible','off','Position',[100 100 1000 300]);
ax  = axes(fig);
imagesc(ax, dateNumeric_est, 1:cfg.horizons, squeeze(pred_q(:, qt_idx(4), :))');
colorbar(ax); colormap(ax,'parula');
datetick(ax, 'x', 'yyyy', 'keepticks');
title(ax, 'BVAR Inflation — Median (p50) by origin & horizon');
xlabel(ax, 'Forecast origin'); ylabel(ax, 'Horizon h (months)');
saveFig(fig, figDir, 'BVAR_infl_median_heatmap.png');

fprintf('\nAll done. Results in: %s\n', outDir);

%% ════════════════════════════════════════════════════════════════════════
%%  LOCAL HELPER
%% ════════════════════════════════════════════════════════════════════════

function plotFanBands(ax, dates, Q)
%PLOTFANBANDS  Draw shaded quantile fan bands.
%   dates : 1×T or T×1 datetime vector
%   Q     : T×7 matrix  [Q05 Q10 Q25 Q50 Q75 Q90 Q95]
    d = dates(:)';
    bands = { [1,7], [0.85 0.9 1],  '5^{th}–95^{th}'; ...
              [2,6], [0.65 0.8 1],  '10^{th}–90^{th}'; ...
              [3,5], [0.4  0.6 1],  '25^{th}–75^{th}' };
    for b = 1:3
        lo = Q(:, bands{b,1}(1))';  hi = Q(:, bands{b,1}(2))';
        xpoly = [d, fliplr(d)];  ypoly = [lo, fliplr(hi)];
        if all(isnan(ypoly)), continue; end
        fill(ax, xpoly, ypoly, bands{b,2}, ...
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

function mkdirs(dirs)
%MKDIRS  Create directories silently if they do not exist.
    for k = 1:numel(dirs)
        if ~exist(dirs{k}, 'dir')
            mkdir(dirs{k});
        end
    end
end
