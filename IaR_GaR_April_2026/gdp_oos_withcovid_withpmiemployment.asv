%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GDP — OUT-OF-SAMPLE (RECURSIVE) — QREG LP + SKEW-T (QUARTERLY)
%  - FAST VERSION: skew-t fitted only at last origin
%  - WITH COVID dummy (same logic as in-sample GDP code)
%  - WITHOUT pmi_composite (4 predictors + constant + covid dummy)
%  - Horizons = 13 (h=1 "current", h=2..13 = 1..12 quarters ahead)
%  - Outputs:
%       * qreg coeffs + last-vintage bootstrap draws
%       * predicted quantiles (all OOS origins)
%       * skew-t params + moments (LAST ORIGIN ONLY)
%       * rolling fan charts (aligned to target dates)
%       * last-origin forward fans (1-year and 3-year)
%       * historical decomposition
%       * Lopez-Salido style Figure 1
%       * skew-t parameter export for sharing
%
% This version: fast + covid / 2026
% This version: 15/04/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS AND PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultAxesFontName','Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1);
rng('default');

baseDir = '\\ma\data\Cross Divisional Work\RASS\IaR Results\IaR_GaR_April_2026\';
cd(baseDir);

addpath(fullfile(baseDir,'intermediate_codes'));
addpath(fullfile(baseDir,'functions'));
addpath(fullfile(baseDir,'Codes','functions','azzalini'));
addpath(fullfile(baseDir,'Codes','functions','CRPS'));
addpath(fullfile(baseDir,'Codes','functions','Simon_qreg'));
addpath(fullfile(baseDir,'Codes','functions','two_piece_normal'));

outputFolder  = fullfile(baseDir,'Outputs');
dropboxFolder = fullfile(baseDir,'DropBoxSink');

ensure = @(p) (exist(p,'dir') || mkdir(p));
ensure(outputFolder);
ensure(fullfile(outputFolder,'sktparam'));
ensure(fullfile(outputFolder,'econ_interpretation_charts'));
ensure(fullfile(dropboxFolder,'predictive_densities'));

% Data file
% fullFileName = fullfile(baseDir,'GaRDataRaw_quarterly_BIS_February_balanced.xlsx');
fullFileName = fullfile(baseDir,'GaRDataRaw_quarterly_BIS_march.xlsx');


%%%%%%%%%%%%%%%%%%
%% CONTROL PANEL
%%%%%%%%%%%%%%%%%%

% startT = datenum(1989,12,31);
startT = datenum(1980,06,30);
% endT   = datenum(2025,12,31);
endT   = datenum(2026,03,31);

ctrynames = {'UK'};
onlyuk    = 1;

% Variables (quarterly GDP) — NO pmi_composite
dep_var_name       = {'g4rgdp'};
current_act_cat    = {'mgdp_yoy'};
leverage_cat       = {'global_credit'};
fci_cat            = {'ciss_uk'};
macro_cond_cat     = {'g4_import_deflator_fuel'};
% lab_mkt_cat        = {'pmi_employment'};
lab_mkt_cat        = {'labour'};


% Skew-t only
model_selection = 2;
modellist = {'ols','skewt','semi-param','two-piece-normal'};

% LP / Quantiles
horizons       = 13;
quantilelevels = 0.05:0.05:0.95;

% First forecast origin
StartEst = datenum(2004,03,31);
% StartEst = datenum(2014,03,31);

% COVID dummy
covid        = 1;
dummyvarname = {'covid'};
covid_dates  = datenum(2020,06,30);   % quarterly COVID date

% GDP option
h_step_gdp = 0;   % 0 = yoy

% Bootstrap (only last origin)
bstOptions.blocksize = 8;
bstOptions.nboot     = 10;    % 5000 in production
bstOptions.ci        = 68;

% Plot horizons
h_list_plot = [2, 5, 9];

%%%%%%%%%%%%%%%%%%%%%%
%% MODEL COMBINATIONS
%%%%%%%%%%%%%%%%%%%%%%
% 
% [i0,i1,i2,i3,i4] = ndgrid(1:size(dep_var_name,2), ...
%                            1:size(current_act_cat,2), ...
%                            1:size(leverage_cat,2), ...
%                            1:size(fci_cat,2), ...
%                            1:size(macro_cond_cat,2));
% 
% vars = {dep_var_name(i0(:)), current_act_cat(i1(:)), leverage_cat(i2(:)), ...
%         fci_cat(i3(:)), macro_cond_cat(i4(:))};

[i0,i1,i2,i3,i4,i5] = ndgrid(1:size(dep_var_name,2), ...
                               1:size(current_act_cat,2), ...
                               1:size(leverage_cat,2), ...
                               1:size(fci_cat,2), ...
                               1:size(macro_cond_cat,2), ...
                               1:size(lab_mkt_cat,2));

vars = {dep_var_name(i0(:)), current_act_cat(i1(:)), leverage_cat(i2(:)), ...
        fci_cat(i3(:)), macro_cond_cat(i4(:)), lab_mkt_cat(i5(:))};


for k = 1:numel(vars)
    v = vars{k};
    if size(v,2) > 1, vars{k} = v.'; end
end

% combo_specifications = [vars{1}, vars{2}, vars{3}, vars{4}, vars{5}];
combo_specifications = [vars{1}, vars{2}, vars{3}, vars{4}, vars{5}, vars{6}];

nSpec = size(combo_specifications,1);

%%%%%%%%%%%%%%%%%%%%%%
%% CONTAINERS
%%%%%%%%%%%%%%%%%%%%%%
coeffqr_OOS             = [];
bootstrapqrg_OOS        = [];
predicted_quantiles_OOS = [];
quarters_origin         = [];
idx_estimation          = [];
idx_end                 = [];
dateNumeric_full        = [];
explvar_all_spec        = [];
exovar_save             = [];

tic;

%%%%%%%%%%%%%%%%%%%%%%
%% LOOP ACROSS SPECS
%%%%%%%%%%%%%%%%%%%%%%

for spec = 1:nSpec

    fprintf('\n=== SPEC %d/%d ===\n', spec, nSpec);

    % Variables expected by dataloading_quarterly
    varnames     = combo_specifications(spec,:);
    lagctryvar   = ones(size(varnames));
    lagglobalvar = [];
    lags         = [lagctryvar lagglobalvar];

    % Load data
    dataloading_quarterly;

    dateNumeric_full = dateNumeric;

    idx_estimation = find(dateNumeric == StartEst, 1, 'first');
    if isempty(idx_estimation)
        error('StartEst (%s) not found in dateNumeric.', datestr(StartEst));
    end
    idx_end     = numel(dateNumeric);
    last_origin = idx_end;

    nOrigins = last_origin - idx_estimation + 1;
    Qn       = numel(quantilelevels);

    % Lagged variable names (BoE info-set convention)
    varnames_lag = strcat('l1', combo_specifications(spec,:));
    depvarname   = varnames_lag(1,1);
    explvarname  = varnames_lag(1,2:end);

    if onlyuk ~= 1, error('UK-only. Set onlyuk=1.'); end

    Ttab    = countryData.UK;
    depvar  = Ttab.(depvarname{1});
    explvar = table2array(Ttab(:, explvarname));

    % COVID dummy: load from countryData and lag by 1 (BoE info set)
    exovar_raw = Ttab.(dummyvarname{1});             % T x 1 (0/1)
    exovar     = [0; exovar_raw(1:end-1)];           % lag by one period

    % Save for later decomposition
    explvar_all_spec(:,:,spec) = explvar;
    exovar_save = exovar;

    % Allocate once
    if isempty(coeffqr_OOS)
        K = size(explvar,2);
        % V = constant + K predictors + covid dummy
        V = 1 + K + 1;

        coeffqr_OOS             = NaN(V, Qn, horizons, nOrigins, nSpec);
        bootstrapqrg_OOS        = NaN(V, Qn, horizons, bstOptions.nboot, nSpec);
        predicted_quantiles_OOS = NaN(nOrigins, Qn, horizons, nSpec);

        StartEstDT   = datetime(StartEst,'ConvertFrom','datenum');
        lastOriginDT = datetime(dateNumeric(last_origin),'ConvertFrom','datenum');
        quarters_origin = StartEstDT:calquarters(1):lastOriginDT;

        dateNumeric_full = dateNumeric;
    end

    % Index of the COVID date in the full sample
    idx_covid = find(dateNumeric == min(covid_dates), 1, 'first');

    % -----------------------------------------------------------------
    % RECURSIVE OOS LOOP
    % -----------------------------------------------------------------
    for endtime = idx_estimation:last_origin

        tOOS = endtime - idx_estimation + 1;
        fprintf('  origin %d/%d  (%s)\n', tOOS, nOrigins, datestr(dateNumeric(endtime)));

        X = explvar(1:endtime,:);
        y = depvar(1:endtime,:);

        % LP targets
        Y_LP = NaN(size(X,1), 1, horizons);
        if h_step_gdp == 0
            for h = 1:horizons
                if size(X,1) > h
                    Y_LP(1:size(X,1)-h, 1, h) = y(1+h:end);
                end
            end
        else
            for h = 1:horizons
                if size(y,1) > h
                    tmp = 100*(log(y(1+h:end)) - log(y(1:end-h)));
                    Y_LP(1:numel(tmp), 1, h) = tmp * (4/h);
                end
            end
        end

        % COVID logic: include dummy only if estimation window covers COVID
        if isempty(idx_covid) || endtime < idx_covid
            covid_qreg = 0;
            exo_now    = 0;
            Xmat       = [ones(size(X,1),1), X];
        else
            covid_qreg = 1;
            exo_now    = exovar(1:endtime,:);
            Xmat       = [ones(size(X,1),1), X];
            % Note: exo_now is passed separately to the function
        end

        bstNow = (endtime == last_origin);

        [bQR, bQRbst] = qfe_qr_local_projection_SL_final( ...
            Y_LP, Xmat, ...
            quantilelevels, (1:horizons)', bstNow, bstOptions, covid_qreg, exo_now);

        % bQR size depends on whether COVID was included:
        % without COVID: (1+K) x Q x H
        % with COVID:    (1+K+1) x Q x H
        nCoeff = size(bQR,1);
        coeffqr_OOS(1:nCoeff,:,:,tOOS,spec) = bQR;

        if bstNow
            bootstrapqrg_OOS(1:nCoeff,:,:,:,spec) = bQRbst;
        end

        % Predict quantiles at this origin
        % Build x_now to match the coefficient vector
        if covid_qreg == 0
            x_now = [1, X(endtime,:)];
        else
            x_now = [1, X(endtime,:), exovar(endtime)];
        end
        % Pad if needed (x_now must be 1 x nCoeff)
        x_now_padded = zeros(1, V);
        x_now_padded(1:nCoeff) = x_now;

        for h = 1:horizons
            predicted_quantiles_OOS(tOOS,:,h,spec) = x_now_padded(1:nCoeff) * bQR(:,:,h);
        end
    end

    % Monotonicity
    predicted_quantiles_OOS(:,:,:,spec) = sort(predicted_quantiles_OOS(:,:,:,spec), 2);

    % -----------------------------------------------------------------
    % FIT SKEW-T — LAST ORIGIN ONLY (fast)
    % -----------------------------------------------------------------
    t_last = nOrigins;
    fprintf('  skew-t fit: last origin only (t=%d)\n', t_last);

    PQ = reshape(predicted_quantiles_OOS(t_last,:,:,spec), 1, 1, Qn, horizons, 1);

    [aa, ~, ~, ~, bb] = fit_skewt_to_quantiles_all( ...
        PQ, quantilelevels, [0.05 0.25 0.5 0.75 0.95], 1, 1:horizons);

    param_skt_last = squeeze(aa);            % 4 x H
    fcstmean_last  = squeeze(bb(:,1,:,:))';  % 1 x H
    fcststdev_last = sqrt(squeeze(bb(:,2,:,:))');
    fcstskew_last  = squeeze(bb(:,3,:,:))';

    % Save actual GDP + meta
    actual_var_long = countryData.(ctrynames{1}).(dep_var_name{1});
    save(fullfile(outputFolder,'actual_gdp_yoy_OOS.mat'), ...
        'actual_var_long','dateNumeric_full','idx_estimation','idx_end', ...
        'StartEst','endT','startT','last_origin');

end % spec loop

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE CORE OBJECTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(fullfile(outputFolder,'qreg_results_gdp_OOS.mat'),           'coeffqr_OOS');
save(fullfile(outputFolder,'predicted_quantiles_gdp_OOS.mat'),     'predicted_quantiles_OOS');
save(fullfile(outputFolder,'bootstrap_results_gdp_OOS.mat'),       'bootstrapqrg_OOS');
save(fullfile(outputFolder,'explanatoryvar_gdp_OOS.mat'),          'explvar_all_spec','exovar_save');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROLLING FAN CHARTS — ALIGNED TO TARGET DATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_to_use = 1;

actualDT = datetime(dateNumeric_full,'ConvertFrom','datenum');

quantilesplot = [0.05 0.10 0.25 0.50 0.75 0.90 0.95];
idx_qt = arrayfun(@(q) find(abs(quantilelevels-q)<1e-8,1,'first'), quantilesplot);
bands  = {'5^{th}–95^{th}','10^{th}–90^{th}','25^{th}–75^{th}'};

T_plot = size(predicted_quantiles_OOS,1);

for ii = 1:numel(h_list_plot)

    h_plot = h_list_plot(ii);

    origin_dates = quarters_origin(1:T_plot);
    fcst_dates   = origin_dates + calquarters(h_plot - 1);

    Q05 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(1), h_plot, spec_to_use));
    Q10 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(2), h_plot, spec_to_use));
    Q25 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(3), h_plot, spec_to_use));
    Q50 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(4), h_plot, spec_to_use));
    Q75 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(5), h_plot, spec_to_use));
    Q90 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(6), h_plot, spec_to_use));
    Q95 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(7), h_plot, spec_to_use));

    % Actual matched to target dates
    [tf, loc] = ismember(datenum(fcst_dates), dateNumeric_full);
    act_plot = NaN(size(fcst_dates));
    act_plot(tf) = actual_var_long(loc(tf));

    figure('Units','normalized','Position',[0.2 0.15 0.6 0.7],'Color','w');
    ax = axes; hold(ax,'on');

    fill([fcst_dates, fliplr(fcst_dates)], [Q05', fliplr(Q95')], [0.85 0.9 1], ...
         'EdgeColor','none','FaceAlpha',1,'DisplayName',bands{1});
    fill([fcst_dates, fliplr(fcst_dates)], [Q10', fliplr(Q90')], [0.65 0.8 1], ...
         'EdgeColor','none','FaceAlpha',1,'DisplayName',bands{2});
    fill([fcst_dates, fliplr(fcst_dates)], [Q25', fliplr(Q75')], [0.4 0.6 1], ...
         'EdgeColor','none','FaceAlpha',1,'DisplayName',bands{3});

    plot(ax, fcst_dates, Q50, '--', 'Color',[0 0 0.7], 'LineWidth',1.2,'DisplayName','50^{th}');
    plot(ax, fcst_dates, act_plot, 'k', 'LineWidth',1.25,'DisplayName','Outturn');

    yline(ax,0,'k-','LineWidth',0.75,'HandleVisibility','off');
    grid(ax,'on');

    ticks = datetime(year(fcst_dates(1)),1,1):calyears(2):fcst_dates(end);
    ax.XTick = ticks; ax.XAxis.TickLabelFormat = 'yyyy';

    legend(ax,'show','Location','northoutside','Orientation','horizontal'); legend boxoff;
    set(ax,'FontSize',14);
    title(ax, sprintf('GDP YoY Forecast (OOS) — h = %d (%dQ ahead)', h_plot, h_plot-1));
    hold(ax,'off');

    print(gcf, fullfile(dropboxFolder,'predictive_densities', ...
        sprintf('gdp_fanchart_OOS_h%d.png', h_plot)), '-dpng', '-r200');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAST-ORIGIN FORWARD FAN: 1–12 QUARTERS AHEAD (3 years)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_to_use = 1;

t_last           = size(predicted_quantiles_OOS, 1);
last_origin_date = quarters_origin(t_last);
fcst_dates_fwd   = last_origin_date + calquarters(1:12);

Q05 = NaN(1,12); Q10 = Q05; Q25 = Q05; Q50 = Q05; Q75 = Q05; Q90 = Q05; Q95 = Q05;

for h_fwd = 1:12
    h_store = h_fwd + 1;
    Q05(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(1), h_store, spec_to_use);
    Q10(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(2), h_store, spec_to_use);
    Q25(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(3), h_store, spec_to_use);
    Q50(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(4), h_store, spec_to_use);
    Q75(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(5), h_store, spec_to_use);
    Q90(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(6), h_store, spec_to_use);
    Q95(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(7), h_store, spec_to_use);
end

act_series = actual_var_long(:);

figure('Units','normalized','Position',[0.18 0.12 0.64 0.72],'Color','w');
ax = axes; hold(ax,'on');

plot(ax, actualDT, act_series, 'k', 'LineWidth', 1.25, 'DisplayName', 'Outturn (GDP YoY)');
xline(ax, last_origin_date, 'k--', 'LineWidth', 0.9, 'DisplayName', 'Last origin');

fill([fcst_dates_fwd, fliplr(fcst_dates_fwd)], [Q05, fliplr(Q95)], [0.85 0.9 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','5^{th}–95^{th}');
fill([fcst_dates_fwd, fliplr(fcst_dates_fwd)], [Q10, fliplr(Q90)], [0.65 0.8 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','10^{th}–90^{th}');
fill([fcst_dates_fwd, fliplr(fcst_dates_fwd)], [Q25, fliplr(Q75)], [0.4 0.6 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','25^{th}–75^{th}');

plot(ax, fcst_dates_fwd, Q50, '--', 'Color', [0 0 0.7], 'LineWidth', 1.2, 'DisplayName','Median (50^{th})');
yline(ax, 0, 'k-', 'LineWidth', 0.75, 'HandleVisibility','off');

grid(ax, 'on');
x_start = datetime(2015,1,1);
xlim(ax, [x_start, fcst_dates_fwd(end)+calquarters(1)]);
ticks = x_start:calyears(2):fcst_dates_fwd(end)+calyears(1);
ax.XTick = ticks; ax.XAxis.TickLabelFormat = 'yyyy';
ylim(ax, 'auto');

legend(ax,'show','Location','northoutside','Orientation','horizontal'); legend boxoff;
set(ax,'FontSize',13);
title(ax, sprintf('GDP YoY — last-origin fan (%s): 1–12 quarters ahead', ...
    datestr(last_origin_date,'QQ yyyy')));
hold(ax,'off');

print(gcf, fullfile(dropboxFolder,'predictive_densities', ...
    'gdp_fanchart_OOS_lastOrigin_12Q.png'), '-dpng', '-r200');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAST-ORIGIN FORWARD FAN: 1 YEAR AHEAD (4 QUARTERS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_to_use = 1;

fcst_dates_1y = last_origin_date + calquarters(1:4);

Q05_1y = NaN(1,4); Q10_1y = Q05_1y; Q25_1y = Q05_1y; Q50_1y = Q05_1y;
Q75_1y = Q05_1y;   Q90_1y = Q05_1y; Q95_1y = Q05_1y;

for h_fwd = 1:4
    h_store = h_fwd + 1;
    Q05_1y(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(1), h_store, spec_to_use);
    Q10_1y(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(2), h_store, spec_to_use);
    Q25_1y(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(3), h_store, spec_to_use);
    Q50_1y(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(4), h_store, spec_to_use);
    Q75_1y(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(5), h_store, spec_to_use);
    Q90_1y(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(6), h_store, spec_to_use);
    Q95_1y(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(7), h_store, spec_to_use);
end

figure('Units','normalized','Position',[0.18 0.12 0.64 0.72],'Color','w');
ax = axes; hold(ax,'on');

plot(ax, actualDT, act_series, 'k', 'LineWidth', 1.25, 'DisplayName', 'Outturn (GDP YoY)');
xline(ax, last_origin_date, 'k--', 'LineWidth', 0.9, 'DisplayName', 'Last origin');

fill([fcst_dates_1y, fliplr(fcst_dates_1y)], [Q05_1y, fliplr(Q95_1y)], [0.85 0.9 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','5^{th}–95^{th}');
fill([fcst_dates_1y, fliplr(fcst_dates_1y)], [Q10_1y, fliplr(Q90_1y)], [0.65 0.8 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','10^{th}–90^{th}');
fill([fcst_dates_1y, fliplr(fcst_dates_1y)], [Q25_1y, fliplr(Q75_1y)], [0.4 0.6 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','25^{th}–75^{th}');

plot(ax, fcst_dates_1y, Q50_1y, '--', 'Color', [0 0 0.7], 'LineWidth', 1.2, 'DisplayName','Median (50^{th})');
yline(ax, 0, 'k-', 'LineWidth', 0.75, 'HandleVisibility','off');

grid(ax, 'on');
x_start = datetime(2015,1,1);
xlim(ax, [x_start, fcst_dates_1y(end)+calquarters(1)]);
ticks = x_start:calyears(2):fcst_dates_1y(end)+calyears(1);
ax.XTick = ticks; ax.XAxis.TickLabelFormat = 'yyyy';
ylim(ax, 'auto');

legend(ax,'show','Location','northoutside','Orientation','horizontal'); legend boxoff;
set(ax,'FontSize',13);
title(ax, sprintf('GDP YoY — last-origin fan (%s): 1–4 quarters ahead (1 year)', ...
    datestr(last_origin_date,'QQ yyyy')));
hold(ax,'off');

print(gcf, fullfile(dropboxFolder,'predictive_densities', ...
    'gdp_fanchart_OOS_lastOrigin_1Y.png'), '-dpng', '-r200');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HISTORICAL DECOMPOSITION (OOS)
%  Includes COVID dummy column in the predictor matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_to_use = 1;

Xspec_full = explvar_all_spec(:,:,spec_to_use);
Xorig      = Xspec_full(idx_estimation:last_origin, :);   % nOrigins x K
exo_orig   = exovar_save(idx_estimation:last_origin, :);   % nOrigins x 1
nOrigins   = size(Xorig,1);

% Full predictor matrix at each origin: [constant, predictors, covid_dummy]
Xfull = [ones(nOrigins,1), Xorig, exo_orig];
V     = size(Xfull,2);

varnames_plot = [{'cons'}, combo_specifications(spec_to_use,2:end), {'covid'}];
varnames_plot = cellfun(@(s) strrep(s,'_',' '), varnames_plot, 'uni',0);

q_to_plot = [0.05, 0.25 0.50 0.75, 0.90];
idx_qt_d  = arrayfun(@(q) find(abs(quantilelevels-q)<1e-8,1), q_to_plot);

hSel  = h_list_plot;
nHsel = numel(hSel);

contribution_pred = NaN(nOrigins, numel(quantilelevels), nHsel, V);
predicted_check   = NaN(nOrigins, numel(quantilelevels), nHsel);

for t = 1:nOrigins
    for ih = 1:nHsel
        h = hSel(ih);
        B_t_h = coeffqr_OOS(1:V,:,h,t,spec_to_use);   % V x Q
%         B_t_h(isnan(B_t_h)) = 0;                        % pre-COVID origins: treat missing COVID coeffs as zero
        predicted_check(t,:,ih) = Xfull(t,:) * B_t_h;
        for v = 1:V
            contribution_pred(t,:,ih,v) = Xfull(t,v) * B_t_h(v,:);
        end
    end
end

quarters_decomp = quarters_origin(1:nOrigins);

for qi = 1:numel(idx_qt_d)
    q_idx = idx_qt_d(qi);

    for ih = 1:nHsel
        h_plot = hSel(ih);

        figure('Units','normalized','Position',[0.2 0.15 0.6 0.7],'Color','w');
        ax = axes; hold(ax,'on');

        C = squeeze(contribution_pred(:, q_idx, ih, :));
        b = bar(ax, quarters_decomp, C, 'stacked', 'BarWidth', 0.7);

        cmap = lines(V);
        for j = 1:V, b(j).FaceColor = cmap(j,:); end

        qline = squeeze(predicted_check(:, q_idx, ih));
        plot(ax, quarters_decomp, qline, 'k-', 'LineWidth',1.3, ...
            'DisplayName', sprintf('q=%.2f', q_to_plot(qi)));

        yline(ax,0,'k-','LineWidth',0.75,'HandleVisibility','off');
        grid(ax,'on');

        ticks = quarters_decomp(1):calyears(2):quarters_decomp(end);
        ax.XTick = ticks; ax.XAxis.TickLabelFormat = 'yyyy';
        set(ax,'FontSize',13);

        lg = legend(ax, [b, findobj(ax,'DisplayName',sprintf('q=%.2f',q_to_plot(qi)))], ...
            [varnames_plot, {sprintf('q=%.2f',q_to_plot(qi))}], ...
            'Location','southoutside','Orientation','horizontal');
        lg.Box = 'off';

        title(ax, sprintf('Historical decomposition (OOS) — h=%d, q=%.0f^{th}', ...
            h_plot, q_to_plot(qi)*100));
        hold(ax,'off');

        print(gcf, fullfile(outputFolder,'econ_interpretation_charts', ...
            sprintf('gdp_decomp_OOS_h%d_q%d.png', h_plot, round(q_to_plot(qi)*100))), ...
            '-dpng', '-r250');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOPEZ-SALIDO FIGURE (LAST VINTAGE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_to_use = 1;

t_last_idx = size(coeffqr_OOS,4);
B_last   = squeeze(coeffqr_OOS(:,:,:,t_last_idx,spec_to_use));   % V x Q x H
std_qreg = squeeze(std(bootstrapqrg_OOS(:,:,:,:,spec_to_use),0,4));

quantilesirf = [0.05 0.1 0.25 0.5 0.75 0.90 0.95];
idx_qt_ls = arrayfun(@(q) find(abs(quantilelevels-q)<1e-8,1,'first'), quantilesirf);
labels_ls = arrayfun(@(q) sprintf('%.0f^{th}',q*100), quantilesirf, 'UniformOutput', false);
colors_ls = lines(numel(idx_qt_ls));

hor_of_interest = [2, 5, 9];
% 4 predictors (no pmi) + covid dummy = 5 rows after constant
% var_order = [1, 2, 3, 4, 5];
var_order = [1, 2, 3, 4, 5, 6];
% 
% varnames_to_plot = { sprintf('ECONOMIC\nACTIVITY'),  sprintf('LEVERAGE\nGROWTH'), ...
%                      sprintf('FINANCIAL\nCONDITIONS'), sprintf('NOMINAL\nINDICATOR'), ...
%                      sprintf('COVID\nDUMMY')};

varnames_to_plot = { sprintf('ECONOMIC\nACTIVITY'),  sprintf('LEVERAGE\nGROWTH'), ...
                     sprintf('FINANCIAL\nCONDITIONS'), sprintf('NOMINAL\nINDICATOR'), ...
                     sprintf('PMI\nEMPLOYMENT'), sprintf('COVID\nDUMMY')};

% fig = figure('Units','centimeters','Position',[1 1 22 24],'Color','w');
fig = figure('Units','centimeters','Position',[1 1 22 28],'Color','w');
tiledlayout(numel(var_order), numel(hor_of_interest), 'TileSpacing','Compact','Padding','Compact');

for k = 1:numel(var_order)
    ii    = var_order(k);
    v_idx = 1 + ii;   % +1 for constant

    for iH = 1:numel(hor_of_interest)
        hor = hor_of_interest(iH);
        ax  = nexttile((k-1)*numel(hor_of_interest) + iH);
        hold(ax,'on'); yline(ax,0,'k');

        for dd = 1:numel(idx_qt_ls)
            q_idx  = idx_qt_ls(dd);
            coeffs = B_last(v_idx, q_idx, hor);
            sds    = std_qreg(v_idx, q_idx, hor);

            errbar = errorbar(ax, dd, coeffs, sds, 's', 'LineWidth',1.5, 'CapSize',0);
            errbar.MarkerFaceColor = colors_ls(dd,:);
            errbar.MarkerEdgeColor = 'none';
        end

        ax.XLim = [0.5, numel(idx_qt_ls)+0.5];
        ax.XTick = 1:numel(idx_qt_ls);
        ax.XTickLabel = labels_ls;
        ax.YAxis.TickLabelFormat = '%.1f';

        if iH == 1, ylabel(ax, varnames_to_plot{ii}, 'Interpreter','none'); end

        if k == 1
            title(ax, sprintf('h=%d (%dQ)', hor, hor-1), ...
                'Units','normalized', 'Position',[0.5,1.10,0]);
        end

        ax.FontSize = 10; ax.Title.FontSize = 12;
        box(ax,'on'); hold(ax,'off');
    end
end

set(fig, 'PaperUnits','centimeters','PaperPosition',[0 0 22 28]);
print(fig, fullfile(outputFolder,'econ_interpretation_charts', ...
    'econ_interpr_gdp_OOS.png'), '-dpng', '-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE LAST-ORIGIN SKEW-T PARAMS (for sharing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_to_use = 1;

% param_skt_last is 4 x H; drop h=1 ("current") -> h=2:13
lc_qtr = param_skt_last(1, 2:end)';
sc_qtr = param_skt_last(2, 2:end)';
sh_qtr = param_skt_last(3, 2:end)';
df_qtr = param_skt_last(4, 2:end)';

fcstmean_qtr  = fcstmean_last(2:end)';
fcststdev_qtr = fcststdev_last(2:end)';
fcstskew_qtr  = fcstskew_last(2:end)';

last_origin_date = quarters_origin(end);
fcst_dates_qtr   = (last_origin_date + calquarters(1:12))';

out_name = fullfile(outputFolder,'sktparam', ...
    sprintf('last_origin_skewt_params_12Q_spec%d.mat', spec_to_use));
save(out_name, ...
    'lc_qtr','sc_qtr','sh_qtr','df_qtr', ...
    'fcstmean_qtr','fcststdev_qtr','fcstskew_qtr', ...
    'fcst_dates_qtr','last_origin_date','spec_to_use');
fprintf('Saved last-origin 12Q skew-t params to:\n  %s\n', out_name);

% Clean export for sharing
outFile = fullfile(outputFolder,'sktparam', ...
    sprintf('quarterly_skewt_params_forSharing_gdp_spec%d.mat', spec_to_use));
save(outFile, 'lc_qtr','sc_qtr','sh_qtr','df_qtr','fcst_dates_qtr', '-v7.3');
fprintf('Saved clean sharing file to:\n  %s\n', outFile);

disp('DONE: Fast GDP OOS pipeline (with COVID, no PMI) complete.');
