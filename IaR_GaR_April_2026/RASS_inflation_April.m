
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INFLATION — OUT-OF-SAMPLE (RECURSIVE) — QUANTILE LP + SKEW-T ONLY (MONTHLY)
%  - Sticks to the OLD script conventions: dataloading_monthly relies on
%    workspace variables (covid, varnames, lags, dummyvarname, fullFileName...)
%  - NO COVID (covid = 0)
%  - ONLY: Quantile LP + Skew-t
%  - Horizons = 13 (1m .. 12m ahead)
%  - Outputs:
%       * qreg coeffs + last-vintage bootstrap draws
%       * predicted quantiles (OOS origins)
%       * skew-t params + moments
%       * fan charts (aligned correctly)
%       * historical decomposition (aligned, for h=[2,4,13] and q=[0.25,0.5,0.75])
%       * Lopez-Salido style Figure 1 (last vintage, h=[2,13])
% 
% This version: 14/01/2026
% This version: 06/02/2026
% This version: 14/04/2026 (This is the full code for oos inflation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS AND PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultAxesFontName','Times');
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1);
rng('default');

baseDir = 'C:\Users\344349\Codes_s\';
cd(baseDir);

addpath(fullfile(baseDir,'Codes'));
addpath(fullfile(baseDir,'Codes','intermediate_codes'));
addpath(fullfile(baseDir,'Codes','functions'));
addpath(fullfile(baseDir,'Codes','functions','azzalini'));
addpath(fullfile(baseDir,'Codes','functions','CRPS'));
addpath(fullfile(baseDir,'Codes','functions','Simon_qreg'));
addpath(fullfile(baseDir,'Codes','functions','two_piece_normal'));  % unused but ok

outputFolder  = fullfile(baseDir,'Outputs');
dropboxFolder = fullfile(baseDir,'DropBoxSink');

ensure = @(p) (exist(p,'dir') || mkdir(p));
ensure(outputFolder);
ensure(fullfile(outputFolder,'sktparam'));
ensure(fullfile(outputFolder,'econ_interpretation_charts'));
ensure(fullfile(dropboxFolder,'predictive_densities'));

% Data file used by dataloading_monthly
%  fullFileName = fullfile(baseDir,'IaRDataRaw_monthly.xlsx'); 
%fullFileName = fullfile(baseDir,'IaRDataRaw_monthly_February.xlsx'); 
 fullFileName = fullfile(baseDir,'IaRDataRaw_monthly_M.xlsx');

%%%%%%%%%%%%%%%%%%
%% CONTROL PANEL
%%%%%%%%%%%%%%%%%%

startT = datenum(1990,02,01);

% endT   = datenum(2025,09,01); gives you results for N25 
% and use endT = datenum(2025,11,01); for F26

% endT   = datenum(2025,09,01);  
% endT   = datenum(2025,11,01);
endT   = datenum(2026,03,01);


 
ctrynames = {'UK'};
onlyuk    = 1;

% Variables (monthly)
dep_var_name           = {'g4cpi'};
infl_persistence_cat   = {'avgcpi'};
infl_expect_cat        = {'infl1_m'};
economic_slack_cat     = {'delta1vu_tightness'};
% supply_shock_cat       = {'yoy_growth_import_defl_FAME'};
supply_shock_cat = {'yoy_growth_import_deflator'};
fci_cat                = {'bond_spread'};




% Skew-t only
model_selection = 2;
modellist = {'ols','skewt','semi-param','two-piece-normal'};

% IMPORTANT: 13 horizons only 
% horizons       = 13;
horizons         = 37;
quantilelevels = 0.05:0.05:0.95;
momentlist     = {'fcstmean','fcststdev','fcstskew'};


% First forecast origin / first estimation vintage
StartEst = datenum(2004,01,01);

% IMPORTANT for the loader: define covid and dummyvarname
covid        = 0;
dummyvarname = {'covid'};   % harmless placeholder (ignored when covid=0)

% Bootstrap: ONLY last origin (set to 5000 for production)
bstOptions.blocksize = 24;
bstOptions.nboot    = 10;    % 5000 in production
bstOptions.ci       = 68;    % placeholder

% Plot horizons
h_list_plot = [2, 4, 13];  % 1m, ~3m, 12m ahead under your h indexing

%%%%%%%%%%%%%%%%%%%%%%
%% MODEL COMBINATIONS
%%%%%%%%%%%%%%%%%%%%%%

[i0,i1,i2,i3,i4,i5] = ndgrid(1:size(dep_var_name,2), ...
                             1:size(infl_persistence_cat,2), ...
                             1:size(infl_expect_cat,2), ...
                             1:size(economic_slack_cat,2), ...
                             1:size(supply_shock_cat,2), ...
                             1:size(fci_cat,2));

vars = {dep_var_name(i0(:)), infl_persistence_cat(i1(:)), infl_expect_cat(i2(:)), ...
        economic_slack_cat(i3(:)), supply_shock_cat(i4(:)), fci_cat(i5(:))};

for k = 1:numel(vars)
    v = vars{k};
    if size(v,2) > 1
        vars{k} = v.';
    end
end

combo_specifications = [vars{1}, vars{2}, vars{3}, vars{4}, vars{5}, vars{6}];
nSpec = size(combo_specifications,1);

%%%%%%%%%%%%%%%%%%%%%%
%% CONTAINERS (allocated after first load)
%%%%%%%%%%%%%%%%%%%%%%
coeffqr_OOS = [];
bootstrapqrg_OOS = [];
predicted_quantiles_OOS = [];
lc_skt = []; sc_skt = []; sh_skt = []; df_skt = [];
fcstmean = []; fcststdev = []; fcstskew = [];
months_origin = [];
idx_estimation = [];
idx_end = [];
dateNumeric_full = [];
explvar_all_spec = [];

tic;

%%%%%%%%%%%%%%%%%%%%%%
%% LOOP ACROSS SPECS
%%%%%%%%%%%%%%%%%%%%%%

for spec = 1:nSpec

    fprintf('\n=== SPEC %d/%d ===\n', spec, nSpec);

    % Variables expected by dataloading_monthly
    varnames = combo_specifications(spec,:);                 % REQUIRED by loader
    lagctryvar   = ones(size(varnames));                     % REQUIRED by loader
    lagglobalvar = [];                                       % REQUIRED by loader
    lags         = [lagctryvar lagglobalvar];                % REQUIRED by loader

    % Load data (creates countryData + dateNumeric, filtered by startT/endT)
    dataloading_monthly;

    % Indices (after loader filtering)
    idx_estimation = find(dateNumeric == StartEst, 1, 'first');
    if isempty(idx_estimation)
        error('StartEst (%s) not found in dateNumeric from dataloading_monthly.', datestr(StartEst));
    end
    idx_end = numel(dateNumeric);

    % Tail guard: keep only origins whose h=2 outturn exists (for aligned fan charts)
% min_h = min(h_list_plot);     % 2
% last_origin = idx_end - min_h; % last origin index in full sample


last_origin = idx_end;


% EXAMPLE (with our data):
% Data runs from 1990m01 to 2025m11.
%
% If h = 13, the forecast made at origin t predicts inflation at t+12.
%
% Example:
%   Origin = 2024m11  →  forecasted date = 2025m11   (OK, observed)
%   Origin = 2024m12  →  forecasted date = 2025m12   (NOT observed)
%
% Therefore, the last valid forecast origin for h = 13 is 2024m11.
% We drop later origins so that all fan charts and outturn lines
% are defined up to the end of the sample.

% max_h = max(h_list_plot);          % longest horizon plotted (e.g. 13)
% last_origin = idx_end - (max_h-1); % last origin with observed y_{t+h}


    if last_origin < idx_estimation
        error('Not enough data between StartEst and endT for the chosen horizons.');
    end

    nOrigins = last_origin - idx_estimation + 1;  % number of forecast origins we will actually produce
    Qn      = numel(quantilelevels);

    % Build lagged names (same as your in-sample)
    varnames_lag = strcat('l1', combo_specifications(spec,:));
    depvarname   = varnames_lag(1,1);
    explvarname  = varnames_lag(1,2:end);

    % UK only
    if onlyuk ~= 1
        error('Panel not implemented in this script. Set onlyuk=1.');
    end

    Ttab   = countryData.UK;
    depvar = Ttab.(depvarname{1});                    % T x 1
    explvar = table2array(Ttab(:, explvarname));      % T x K

    % Save predictors for later decomposition (full sample, we’ll slice later)
    explvar_all_spec(:,:,spec) = explvar; 
    save(fullfile(outputFolder,'explanatoryvar_inflation_OOS.mat'),'explvar_all_spec');

    % Allocate once (first spec only)
    if isempty(coeffqr_OOS)
        K = size(explvar,2);
        V = 1 + K;

        coeffqr_OOS             = NaN(V, Qn, horizons, nOrigins, nSpec);
        bootstrapqrg_OOS        = NaN(V, Qn, horizons, bstOptions.nboot, nSpec); % only last origin filled
        predicted_quantiles_OOS = NaN(nOrigins, Qn, horizons, nSpec);

        lc_skt   = NaN(nOrigins, horizons, nSpec);
        sc_skt   = NaN(nOrigins, horizons, nSpec);
        sh_skt   = NaN(nOrigins, horizons, nSpec);
        df_skt   = NaN(nOrigins, horizons, nSpec);

        fcstmean  = NaN(nOrigins, horizons, nSpec);
        fcststdev = NaN(nOrigins, horizons, nSpec);
        fcstskew  = NaN(nOrigins, horizons, nSpec);

        % Origin date vector (ONLY for produced origins)
        StartEstDT = datetime(StartEst,'ConvertFrom','datenum');
        lastOriginDT = datetime(dateNumeric(last_origin),'ConvertFrom','datenum');
        months_origin = StartEstDT:calmonths(1):lastOriginDT;

        dateNumeric_full = dateNumeric;
    end

    % ---------------------------------------------------------------------
    % RECURSIVE OOS LOOP (origins)
    % ---------------------------------------------------------------------
    for endtime = idx_estimation:last_origin

        tOOS = endtime - idx_estimation + 1;
        fprintf('  origin %d/%d  (%s)\n', tOOS, nOrigins, datestr(dateNumeric(endtime)));

        X = explvar(1:endtime,:);
        y = depvar(1:endtime,:);

        % LP target container (T x 1 x H)
        Y_LP = NaN(size(X,1),1,horizons);
        for h = 1:horizons
            if size(X,1) > h
                Y_LP(1:size(X,1)-h,1,h) = y(1+h:end);
            end
        end

        % Bootstrap only for last origin (huge speed-up)
        bstNow = (endtime == last_origin);

        [bQR, bQRbst] = qfe_qr_local_projection_SL_final( ...
            Y_LP, [ones(size(X,1),1), X], ...
            quantilelevels, (1:horizons)', bstNow, bstOptions, 0, 0);

        % Store coefficients for this origin
        coeffqr_OOS(:,:,:,tOOS,spec) = bQR;

        % Store bootstrap draws only at last origin
        if bstNow
            bootstrapqrg_OOS(:,:,:,:,spec) = bQRbst;
        end

        % Predict quantiles at this origin.  
        % 1 = constant; X(endtime,:) = the predictor values from the last
        % row of the dataset = eg: Nov 2025 
        x_now = [1, X(endtime,:)];
        for h = 1:horizons
            predicted_quantiles_OOS(tOOS,:,h,spec) = x_now * bQR(:,:,h);  % (1xV)*(VxQ); coefficients we estimated are in bQR(:,:,h)
        end
    end

    % Ensure monotonicity across quantiles
    predicted_quantiles_OOS(:,:,:,spec) = sort(predicted_quantiles_OOS(:,:,:,spec), 2);

    % ---------------------------------------------------------------------
    % FIT SKEW-T PARAMS + MOMENTS for produced origins only
    % ---------------------------------------------------------------------
    for t = 1:nOrigins
        for h = 1:horizons
            [lc_skt(t,h,spec), sc_skt(t,h,spec), sh_skt(t,h,spec), df_skt(t,h,spec)] = ...
                QuantilesInterpolation( squeeze(predicted_quantiles_OOS(t,:,h,spec)), quantilelevels );

            [fcstmean(t,h,spec), fcststdev(t,h,spec), fcstskew(t,h,spec), ~] = ...
                moments_skewt_distr_sm_updated(lc_skt(t,h,spec), sc_skt(t,h,spec), ...
                                               sh_skt(t,h,spec), df_skt(t,h,spec));
        end
    end

    % Save actual inflation and meta (full sample from loader)
    actual_var_long = countryData.(ctrynames{1}).(dep_var_name{1});
    save(fullfile(outputFolder,'actual_inflation_mom_OOS.mat'), ...
        'actual_var_long','dateNumeric_full','idx_estimation','idx_end','StartEst','endT','startT','last_origin');

end % spec loop

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE CORE OOS OBJECTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(fullfile(outputFolder,'qreg_results_inflation_OOS.mat'), 'coeffqr_OOS');                     % V×Q×H×T×S
save(fullfile(outputFolder,'predicted_quantiles_inflation_OOS.mat'), 'predicted_quantiles_OOS'); % T×Q×H×S
save(fullfile(outputFolder,'bootstrap_results_inflation_OOS.mat'), 'bootstrapqrg_OOS');          % V×Q×H×B×S (only last origin filled)

save(fullfile(outputFolder,'sktparam','skewtparam_inflation_OOS.mat'), 'lc_skt','sc_skt','sh_skt','df_skt');
save(fullfile(outputFolder,'sktparam','moments_inflation_OOS.mat'),    'fcstmean','fcststdev','fcstskew');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE MOMENTS TO EXCEL (origins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_to_use = 1;

filename = fullfile(outputFolder, sprintf('%s_OUT_OF_SAMPLE_spec_%d_inflation_mom.xlsx', ...
                               modellist{model_selection}, spec_to_use));

for i = 1:numel(momentlist)
    currentVar = eval(momentlist{i});           % T×H×S
    currentVar_final = currentVar(:,:,spec_to_use);

    Ttbl = [table(months_origin', 'VariableNames', {'Dates'}), array2table(currentVar_final)];
    forecastNames = strcat("h_", string(0:horizons-1));
    Ttbl.Properties.VariableNames(2:end) = cellstr(forecastNames);

    writetable(Ttbl, filename, 'Sheet', momentlist{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FAN CHARTS (OOS) — CORRECTLY ALIGNED TO FORECASTED DATES


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPRETATION OF THE OOS FAN CHARTS (IMPORTANT)
%
% In this OOS setup, each fan chart is indexed by FORECASTED (target) dates,
% not by forecast origin dates.
%
% Key conventions in this script:
%
%   - Horizons are indexed as h = 1,...,13
%   - By construction:
%         h = 1  -> "current" (0 months ahead)
%         h = 2  -> 1 month ahead
%         ...
%         h = 13 -> 12 months ahead (≈ 1 year ahead)
%
% The x-axis dates are constructed as:
%
%     origin_dates = months_origin(1:T_plot);
%     fcst_dates   = origin_dates + calmonths(h_plot-1);
%
% Therefore, for a given horizon h_plot:
%
%     fcst_dates(t) = origin_dates(t) + (h_plot-1) months
%
% Each point on the fan represents:
%   - A forecast made at origin date t,
%   - For inflation at the TARGET date t + (h_plot-1) months,
%   - Using only information available up to that origin (recursive OOS).
%
% The shaded bands are the predictive quantile intervals:
%   - 5th–95th percentile band  (widest)
%   - 10th–90th percentile band
%   - 25th–75th percentile band (tightest)
%
% The dashed blue line is the median forecast (50th percentile).
% The black line is the realized outturn, aligned to the SAME target dates.
%
% IMPORTANT ABOUT SAMPLE END:
%
% Let dateNumeric_full end in Nov 2025.
% Origins are restricted by:
%
%     min_h = min(h_list_plot);   % here min_h = 2 (1 month ahead)
%     last_origin = idx_end - min_h;
%
% So the LAST origin is Sep 2025.
%
% For the 1-year-ahead fan (h_plot = 13):
%   - Target date = origin + 12 months
%   - Last target date = Sep 2025 + 12 months = Sep 2026
%
% Hence, the 1-year fan EXTENDS BEYOND THE DATA END DATE.
%
% For target dates after Nov 2025:
%   - No outturn is available,
%   - act_plot remains NaN,
%   - The black "Outturn" line stops,
%   - The fan shows PURE FORECASTS only.
%
% If instead we want ALL fan points to have realized outturns,
% we must restrict origins using the MAX horizon, not the MIN horizon:
%
%     max_h = max(h_list_plot);      % e.g. 13
%     last_origin = idx_end - max_h; % ensures fcst_dates never exceed endT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(outputFolder,'predicted_quantiles_inflation_OOS.mat'),'predicted_quantiles_OOS');
load(fullfile(outputFolder,'actual_inflation_mom_OOS.mat'),'actual_var_long','dateNumeric_full','StartEst','endT','last_origin');

spec_to_use = 1;

actualDT = datetime(dateNumeric_full,'ConvertFrom','datenum');

quantilesplot = [0.05 0.10 0.25 0.50 0.75 0.90 0.95];
idx_qt = arrayfun(@(q) find(abs(quantilelevels-q)<1e-8,1,'first'), quantilesplot);
bands = {'5^{th}–95^{th}','10^{th}–90^{th}','25^{th}–75^{th}'};

T_plot = size(predicted_quantiles_OOS,1);

for ii = 1:numel(h_list_plot)

    h_plot = h_list_plot(ii);

    origin_dates = months_origin(1:T_plot);
    fcst_dates   = origin_dates + calmonths(h_plot-1);

    Q05 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(1), h_plot, spec_to_use));
    Q10 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(2), h_plot, spec_to_use));
    Q25 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(3), h_plot, spec_to_use));
    Q50 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(4), h_plot, spec_to_use));
    Q75 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(5), h_plot, spec_to_use));
    Q90 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(6), h_plot, spec_to_use));
    Q95 = squeeze(predicted_quantiles_OOS(1:T_plot, idx_qt(7), h_plot, spec_to_use));

    % Actual series matched to fcst_dates
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
    title(ax, sprintf('Inflation Forecast (OOS) — h = %d months ahead', h_plot));
    hold(ax,'off');

    print(gcf, fullfile(dropboxFolder,'predictive_densities', ...
        sprintf('inflation_fanchart_OOS_h%d.png',h_plot)), '-dpng', '-r200');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HISTORICAL DECOMPOSITION (OOS)
% Contributions to fitted quantiles at each origin using that origin's betas
% We compute ONLY for h_list_plot to keep it fast.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_to_use = 1;

load(fullfile(outputFolder,'qreg_results_inflation_OOS.mat'),'coeffqr_OOS');
load(fullfile(outputFolder,'explanatoryvar_inflation_OOS.mat'),'explvar_all_spec');

Xspec_full = explvar_all_spec(:,:,spec_to_use);

% Slice to the SAME origin window used in estimation:
% origins correspond to full-sample indices idx_estimation:last_origin,
% but the predictor row used for forecast at origin endtime is explvar(endtime,:).
% So the origin predictor matrix is:
load(fullfile(outputFolder,'actual_inflation_mom_OOS.mat'),'idx_estimation','last_origin');
Xorig = Xspec_full(idx_estimation:last_origin, :);     % nOrigins × K
nOrigins = size(Xorig,1);

Xfull = [ones(nOrigins,1), Xorig];
V = size(Xfull,2);
Qn = numel(quantilelevels);

varnames_plot = [{'cons'}, combo_specifications(spec_to_use,2:end)];
varnames_plot = cellfun(@(s) strrep(s,'_',' '), varnames_plot, 'uni',0);

% Contribution: t × q × hSel × v
hSel = h_list_plot;
nHsel = numel(hSel);

contribution_pred = NaN(nOrigins, Qn, nHsel, V);
predicted_check   = NaN(nOrigins, Qn, nHsel);

for t = 1:nOrigins
    for ih = 1:nHsel
        h = hSel(ih);
        B_t_h = coeffqr_OOS(:,:,h,t,spec_to_use);   % V×Q
        predicted_check(t,:,ih) = Xfull(t,:) * B_t_h;
        for v = 1:V
            contribution_pred(t,:,ih,v) = Xfull(t,v) * (B_t_h(v,:));
        end
    end
end

q_to_plot = [0.25 0.50 0.75];
idx_qt    = arrayfun(@(q) find(abs(quantilelevels-q)<1e-8,1), q_to_plot);

months_decomp = months_origin(1:nOrigins);

for qi = 1:numel(idx_qt)
    q_idx = idx_qt(qi);

    for ih = 1:nHsel
        h_plot = hSel(ih);

        figure('Units','normalized','Position',[0.2 0.15 0.6 0.7],'Color','w');
        ax = axes; hold(ax,'on');

        C = squeeze(contribution_pred(:, q_idx, ih, :)); % nOrigins×V
        b = bar(ax, months_decomp, C, 'stacked', 'BarWidth', 0.7);

        cmap = lines(V);
        for j=1:V, b(j).FaceColor = cmap(j,:); end

        qline = squeeze(predicted_check(:, q_idx, ih));
        plot(ax, months_decomp, qline, 'k-', 'LineWidth',1.3, ...
            'DisplayName', sprintf('q=%.2f', q_to_plot(qi)));

        yline(ax,0,'k-','LineWidth',0.75,'HandleVisibility','off');
        grid(ax,'on');

        ticks = months_decomp(1):calyears(2):months_decomp(end);
        ax.XTick = ticks; ax.XAxis.TickLabelFormat='yyyy';
        set(ax,'FontSize',13);

        lg = legend(ax, [b, findobj(ax,'DisplayName',sprintf('q=%.2f',q_to_plot(qi)))], ...
            [varnames_plot, {sprintf('q=%.2f',q_to_plot(qi))}], ...
            'Location','southoutside','Orientation','horizontal');
        lg.Box='off';

        title(ax, sprintf('Historical decomposition (OOS) — h=%d, q=%.0f^{th}', ...
            h_plot, q_to_plot(qi)*100));

        hold(ax,'off');

        print(gcf, fullfile(outputFolder,'econ_interpretation_charts', ...
            sprintf('inflation_decomp_OOS_h%d_q%d.png', h_plot, round(q_to_plot(qi)*100))), ...
            '-dpng', '-r250');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOPEZ-SALIDO & LORIA STYLE FIGURE 1 (OOS — last origin, 1m & 1y)
% Uses bootstrap SDs from last origin only (as in the framework).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec_to_use = 1;

load(fullfile(outputFolder,'qreg_results_inflation_OOS.mat'),'coeffqr_OOS');
load(fullfile(outputFolder,'bootstrap_results_inflation_OOS.mat'),'bootstrapqrg_OOS');

t_last = size(coeffqr_OOS,4);                           % last produced origin index
B_last = squeeze(coeffqr_OOS(:,:,:,t_last,spec_to_use)); % V×Q×H

% SDs from bootstrap draws (only last origin is filled)
std_qreg = squeeze(std(bootstrapqrg_OOS(:,:,:,:,spec_to_use),0,4)); % V×Q×H

quantilesirf = [0.05 0.1 0.25 0.5 0.75 0.90 0.95];
idx_qt = arrayfun(@(q) find(abs(quantilelevels-q)<1e-8, 1,'first'), quantilesirf);
labels = arrayfun(@(q) sprintf('%.0f^{th}',q*100), quantilesirf, 'UniformOutput', false);
colors = lines(numel(idx_qt));

hor_of_interest = [2, 13];   % 1m and 1y under horizons=13
var_order = [1,2,3,5,4];     % predictor ordering (excluding constant)

varnames_to_plot = { sprintf('INFLATION\nPERSISTENCE'),  sprintf('INFLATION\nEXPECTATIONS'), ...
                     sprintf('ECONOMIC\nTIGHTNESS'), sprintf('EXTERNAL\nCONDITIONS'), ...
                     sprintf('FINANCIAL\nCONDITIONS')};

fig = figure('Units','centimeters','Position',[1 1 22 20],'Color','w');
t = tiledlayout(5,2,'TileSpacing','Compact','Padding','Compact');

for k = 1:numel(var_order)

    ii = var_order(k);     % predictor index 1..5
    v_idx = 1 + ii;        % +1 for constant column

    for iH = 1:numel(hor_of_interest)

        hor = hor_of_interest(iH);
        ax = nexttile((k-1)*numel(hor_of_interest) + iH);

        hold(ax,'on'); yline(ax,0,'k');

        for dd = 1:numel(idx_qt)
            q_idx = idx_qt(dd);
            coeffs = B_last(v_idx, q_idx, hor);
            sds    = std_qreg(v_idx, q_idx, hor);

            errbar = errorbar(ax, dd, coeffs, sds, 's', ...
                         'LineWidth',1.5, 'CapSize',0);
            errbar.MarkerFaceColor = colors(dd,:);
            errbar.MarkerEdgeColor = 'none';
        end

        ax.XLim = [0.5, numel(idx_qt)+0.5];
        ax.XTick = 1:numel(idx_qt);
        ax.XTickLabel = labels;
        ax.YAxis.TickLabelFormat = '%.1f';

        if iH==1
            ylabel(ax, varnames_to_plot{ii}, 'Interpreter','none');
        end

        if k==1
            if hor == 2
                title(ax, '1-MONTH', 'Units','normalized', 'Position',[0.5,1.10,0]);
            else
                title(ax, '1-YEAR', 'Units','normalized', 'Position',[0.5,1.10,0]);
            end
        end

        ax.FontSize = 10;
        ax.Title.FontSize = 12;
        box(ax,'on');
        hold(ax,'off');

    end
end

set(fig, 'PaperUnits','centimeters','PaperPosition',[0 0 22 20]);
print(fig, fullfile(outputFolder,'econ_interpretation_charts', ...
    'econ_interpr_inflation_OOS_1m_1y.png'), '-dpng', '-r300');

disp('DONE: OOS inflation (horizons=13) skew-t pipeline complete.');


% Each of the 12 points (h = 1..12) is a future-month forecast, made at the last origin (Nov 2025), using the quantile regression coefficients estimated up to that date. That s what forms the 12‑month fan extending into Nov 2026.
 
% -------------------------------------------------------------------------
% FIGURE__last-origin fan (h = 1…12): scenario
%
% The fan uses only the information available at the last origin
% (Nov 2025). At that date, the model takes the predictor vector x_t 
% (all explanatory variables observed in Nov 2025) and combines it with the
% horizon‑specific quantile regression coefficients β_{t,h} estimated using
% data up to Nov 2025.
%
% For each forecast horizon h = 1,…,12, corresponding to target months
% Dec 2025 → Nov 2026, the model computes:
%
%       ŷ_{t+h}^{(q)} = x_t * β_{t,h}^{(q)}
%
% producing a full predictive distribution across quantiles
% (5th, 10th, 25th, median, 75th, 90th, 95th).
%
% These 12 predictive distributions form the shaded fan to the right of the
% vertical line at Nov 2025. They represent the model’s view, as of Nov 2025,
% of inflation over the next 12 months.
% -------------------------------------------------------------------------

% 

%% === Single-origin fan (last origin only) across forward horizons ===
% Desired: h_fwd = 1..12 => y_{t+h_fwd} | info at last origin t
% Mapping to stored horizon index: h_store = h_fwd + 1 (since h_store=1 is "current")

spec_to_use = 1;

% Find the last produced origin index
t_last = size(predicted_quantiles_OOS, 1);

% Dates
actualDT = datetime(dateNumeric_full,'ConvertFrom','datenum');
last_origin_date = months_origin(t_last);                 % e.g., 2025-11-01
fcst_dates = last_origin_date + calmonths(1:12);         % Dec-2025 ... Nov-2026

% Pick quantiles for bands
quantilesplot = [0.05 0.10 0.25 0.50 0.75 0.90 0.95];
idx_qt = arrayfun(@(q) find(abs(quantilelevels-q)<1e-8,1,'first'), quantilesplot);

% Extract last-origin quantiles across forward horizons (This line initializes 'empty 1×12 vectors' for each quantile band so they can later be filled with h=1..12 forecasts from the last origin that form the fan chart.
% vectors store future predictive distributions for each of the next 12
% months. later, in the loop, we fill these vectors with forecasts.

Q05 = NaN(1,12); Q10 = Q05; Q25 = Q05; Q50 = Q05; Q75 = Q05; Q90 = Q05; Q95 = Q05;

% t_last = last origin (Nov 2025)
% h = horizon index (2..13 because of the convention used in the code) the result is a set of quantiles (5th, 10th, 25th, …, 95th percentiles)
% for each future month (Dec 2025 → Nov 2026), the code:takes the regression coefficients estimated using data up to Nov 2025. Multiplies them by the predictor values at Nov 2025. Gets a full predictive distribution (skew‑t fitted to quantiles). Plots those quantiles as a fan.

for h_fwd = 1:12
    h_store = h_fwd + 1;  % map to stored horizon index
    Q05(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(1), h_store, spec_to_use);
    Q10(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(2), h_store, spec_to_use);
    Q25(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(3), h_store, spec_to_use);
    Q50(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(4), h_store, spec_to_use);
    Q75(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(5), h_store, spec_to_use);
    Q90(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(6), h_store, spec_to_use);
    Q95(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(7), h_store, spec_to_use);
end

% So the structure becomes:

% Q05(h) = 5th percentile forecast for horizon h
% Q10(h) = 10th percentile
% Q25(h) = 25th percentile
% Q50(h) = median forecast
% Q75(h) = 75th percentile
% Q90(h) = 90th percentile
% Q95(h) = 95th percentile
% 
% These are the fan bands.

% Each element of the vectors (Q05(h), Q10(h), ...) corresponds to:
% 
% The predicted quantile for month = last\_observed\_month + h
% 
% Example (if last observed = Nov 2025):
% 
% | Index h | Calendar month | Meaning                |
% | ------- | -------------- | ---------------------- |
% | 1       | Dec 2025       | Forecast 1-month ahead |
% | 2       | Jan 2026       | 2-month ahead          |
% | 3       | Feb 2026       | 3-month ahead          |
% | …       | …              | …                      |
% | 12      | Nov 2026       | 12-month ahead         |

% So the vectors store future predictive distributions for each of the next 12 months.


% Actual CPI (plot until last observed)
act_series = actual_var_long(:);
% Create a trimmed version up to the last observed date
last_obs_date = actualDT(end);               % should align with endT (Nov 2025)
tmask = actualDT <= last_obs_date;
actDT_trim = actualDT(tmask);
act_trim   = act_series(tmask);

% Plot
figure('Units','normalized','Position',[0.18 0.12 0.64 0.72],'Color','w');
ax = axes; hold(ax,'on');

% Observed CPI (through Nov 2025)
plot(ax, actDT_trim, act_trim, 'k', 'LineWidth', 1.25, 'DisplayName', 'Outturn (CPI)');

% Vertical line at last origin / last observed data point
xline(ax, last_origin_date, 'k--', 'LineWidth', 0.9, 'DisplayName', 'Last origin');

% Fan bands (only last origin, across horizons)
% Build patches across the 12 forecast target dates
fill([fcst_dates, fliplr(fcst_dates)], [Q05, fliplr(Q95)], [0.85 0.9 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','5^{th}–95^{th}');
fill([fcst_dates, fliplr(fcst_dates)], [Q10, fliplr(Q90)], [0.65 0.8 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','10^{th}–90^{th}');
fill([fcst_dates, fliplr(fcst_dates)], [Q25, fliplr(Q75)], [0.4 0.6 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','25^{th}–75^{th}');

% Median
plot(ax, fcst_dates, Q50, '--', 'Color', [0 0 0.7], 'LineWidth', 1.2, 'DisplayName','Median (50^{th})');

% Zero line
yline(ax, 0, 'k-', 'LineWidth', 0.75, 'HandleVisibility','off');

% Cosmetics
grid(ax, 'on');
% Set x-limits to show both the observed history and the 12 forecast months
xlim(ax, [actDT_trim(max(end-120,1)), fcst_dates(end)+calmonths(1)]); % last ~10 years + forecasts
ticks = datetime(year(actDT_trim(1)),1,1):calyears(2):fcst_dates(end);
ax.XTick = ticks; ax.XAxis.TickLabelFormat = 'yyyy';

legend(ax,'show','Location','northoutside','Orientation','horizontal'); legend boxoff;
set(ax,'FontSize',13);
 title(ax, sprintf('Inflation (CPI) — last-origin fan (Nov 2025 origin): 1–12 months ahead'));
% last_origin_date = months_origin(end);
% title(ax, sprintf('Inflation (CPI) — last-origin fan (%s origin): 1–12 months ahead', ...
%                   upper(datestr(last_origin_date,'mmm yyyy'))));
hold(ax,'off');


%% ===CODE SNIPPET that extracts and saves last-origin (Nov 2025) 12-month-ahead skew-t parameters ===
% We take the last produced origin (t_last) and horizons h=2:13 (1..12 months ahead)
% and save: lc_last12, sc_last12, sh_last12, df_last12 + their target dates.
% Horizon mapping: Because the arrays use h=1 as "current", the 12 forward months are h=2:13.

clearvars -except baseDir outputFolder dropboxFolder months_origin quantilelevels

spec_to_use = 1;  % adjust if needed

% Load the parameter and meta files
load(fullfile(outputFolder,'sktparam','skewtparam_inflation_OOS.mat'), 'lc_skt','sc_skt','sh_skt','df_skt');
load(fullfile(outputFolder,'actual_inflation_mom_OOS.mat'), 'StartEst','endT','last_origin'); 
% 'months_origin' was created in-memory in your main script. If not in workspace:
if ~exist('months_origin','var')
    % Rebuild months_origin from meta if necessary
    load(fullfile(outputFolder,'predicted_quantiles_inflation_OOS.mat'),'predicted_quantiles_OOS');
    % Recover StartEst from meta and infer months_origin from length T
    T_plot = size(predicted_quantiles_OOS,1);
    StartEstDT = datetime(StartEst,'ConvertFrom','datenum');
    % last_origin date is available via 'last_origin' + 'dateNumeric_full' if needed
    % but we only need StartEstDT and T_plot to rebuild months_origin:
    months_origin = StartEstDT + calmonths(0:T_plot-1);
end

% Determine indices
t_last = size(lc_skt, 1);   % last produced origin (Nov 2025)
h_idx  = 2:13;              % 1..12 months ahead under your convention

% Handle spec dimension (arrays are [T x H x S] or [T x H] if S=1)
get_slice = @(A) squeeze(A(t_last, h_idx, min(spec_to_use, max(1, size(A,3)))));

if ndims(lc_skt) == 3
    lc_last12 = get_slice(lc_skt);   % 1x12 -> 1x12 or 12x1 after squeeze
    sc_last12 = get_slice(sc_skt);
    sh_last12 = get_slice(sh_skt);
    df_last12 = get_slice(df_skt);
else
    % If no spec dimension (T x H)
    lc_last12 = squeeze(lc_skt(t_last, h_idx));
    sc_last12 = squeeze(sc_skt(t_last, h_idx));
    sh_last12 = squeeze(sh_skt(t_last, h_idx));
    df_last12 = squeeze(df_skt(t_last, h_idx));
end

% Ensure column vectors (nice for saving/inspection)
lc_last12 = lc_last12(:);
sc_last12 = sc_last12(:);
sh_last12 = sh_last12(:);
df_last12 = df_last12(:);

% Build the corresponding target dates: last origin + 1..12 months
last_origin_date = months_origin(t_last);          % e.g., 2025-11-01
fcst_dates       = last_origin_date + calmonths(1:12);
fcst_dates       = fcst_dates(:);

% (Optional) also extract last-origin moments (mean, stdev, skew) for h=2:13
moments_path = fullfile(outputFolder,'sktparam','moments_inflation_OOS.mat');
if exist(moments_path, 'file')
    load(moments_path, 'fcstmean','fcststdev','fcstskew');  % [T x H x S]
    if ndims(fcstmean) == 3
        fcstmean_last12  = squeeze(fcstmean(t_last,  h_idx, spec_to_use));  fcstmean_last12  = fcstmean_last12(:);
        fcststdev_last12 = squeeze(fcststdev(t_last, h_idx, spec_to_use)); fcststdev_last12 = fcststdev_last12(:);
        fcstskew_last12  = squeeze(fcstskew(t_last,  h_idx, spec_to_use));  fcstskew_last12  = fcstskew_last12(:);
    else
        fcstmean_last12  = squeeze(fcstmean(t_last,  h_idx));  fcstmean_last12  = fcstmean_last12(:);
        fcststdev_last12 = squeeze(fcststdev(t_last, h_idx)); fcststdev_last12 = fcststdev_last12(:);
        fcstskew_last12  = squeeze(fcstskew(t_last,  h_idx));  fcstskew_last12  = fcstskew_last12(:);
    end
else
    fcstmean_last12  = [];
    fcststdev_last12 = [];
    fcstskew_last12  = [];
end

% Save all to a single .mat for convenience
out_name = fullfile(outputFolder,'sktparam', sprintf('last_origin_skewt_params_12mo_spec%d.mat', spec_to_use));
save(out_name, ...
    'lc_last12','sc_last12','sh_last12','df_last12', ...
    'fcstmean_last12','fcststdev_last12','fcstskew_last12', ...
    'fcst_dates','last_origin_date','spec_to_use');

fprintf('Saved last-origin 12-month skew-t params to:\n  %s\n', out_name);



%% === Last-origin fan (36 months ahead): Dec-2025 → Nov-2028 ===
% Requires that you have previously saved predicted_quantiles_OOS and have horizons >= 37.
% Mapping: forward horizon h_fwd = 1..36  → stored index h_store = h_fwd + 1 (skip "current").

spec_to_use = 1;

% Load what we need
load(fullfile(outputFolder,'predicted_quantiles_inflation_OOS.mat'),'predicted_quantiles_OOS');
load(fullfile(outputFolder,'actual_inflation_mom_OOS.mat'),'actual_var_long','dateNumeric_full','StartEst'); 

% Rebuild months_origin if not in workspace
if ~exist('months_origin','var')
    T_plot = size(predicted_quantiles_OOS,1);
    StartEstDT = datetime(StartEst,'ConvertFrom','datenum');
    months_origin = StartEstDT + calmonths(0:T_plot-1);
end

% Check that we have enough horizons (need at least 37: 0..36)
H_avail = size(predicted_quantiles_OOS, 3);
if H_avail < 37
    error('Not enough horizons in predicted_quantiles_OOS (have %d, need >= 37). Re-run estimation with horizons = 37.', H_avail);
end

t_last = size(predicted_quantiles_OOS, 1);  % last produced origin (Nov 2025)
last_origin_date = months_origin(t_last);

% Target dates for 36 months ahead
fcst_dates = last_origin_date + calmonths(1:36);

% Quantile indices to plot
quantilesplot = [0.05 0.10 0.25 0.50 0.75 0.90 0.95];
if ~exist('quantilelevels','var')
    % If not in workspace, reconstruct the default used in your script
    quantilelevels = 0.05:0.05:0.95;
end
idx_qt = arrayfun(@(q) find(abs(quantilelevels-q)<1e-8,1,'first'), quantilesplot);

% Extract last-origin quantiles across h_fwd = 1..36 (map to h_store = h_fwd+1)
Q05 = NaN(1,36); Q10 = Q05; Q25 = Q05; Q50 = Q05; Q75 = Q05; Q90 = Q05; Q95 = Q05;
for h_fwd = 1:36
    h_store = h_fwd + 1;  % because h_store=1 is "current"
    Q05(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(1), h_store, spec_to_use);
    Q10(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(2), h_store, spec_to_use);
    Q25(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(3), h_store, spec_to_use);
    Q50(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(4), h_store, spec_to_use);
    Q75(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(5), h_store, spec_to_use);
    Q90(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(6), h_store, spec_to_use);
    Q95(h_fwd) = predicted_quantiles_OOS(t_last, idx_qt(7), h_store, spec_to_use);
end

% Observed CPI (stops at last actual date)
actualDT = datetime(dateNumeric_full,'ConvertFrom','datenum');
act_series = actual_var_long(:);
xline_date = last_origin_date;

% Plot
figure('Units','normalized','Position',[0.18 0.12 0.64 0.72],'Color','w');
ax = axes; hold(ax,'on');

% Observed CPI
plot(ax, actualDT, act_series, 'k', 'LineWidth', 1.25, 'DisplayName', 'Outturn (CPI)');

% Vertical line at last origin
xline(ax, xline_date, 'k--', 'LineWidth', 0.9, 'DisplayName', 'Last origin');

% Fan bands
fill([fcst_dates, fliplr(fcst_dates)], [Q05, fliplr(Q95)], [0.85 0.9 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','5^{th}–95^{th}');
fill([fcst_dates, fliplr(fcst_dates)], [Q10, fliplr(Q90)], [0.65 0.8 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','10^{th}–90^{th}');
fill([fcst_dates, fliplr(fcst_dates)], [Q25, fliplr(Q75)], [0.4 0.6 1], ...
     'EdgeColor','none','FaceAlpha',1, 'DisplayName','25^{th}–75^{th}');

% Median
plot(ax, fcst_dates, Q50, '--', 'Color', [0 0 0.7], 'LineWidth', 1.2, 'DisplayName','Median (50^{th})');

% Zero line and grid
yline(ax, 0, 'k-', 'LineWidth', 0.75, 'HandleVisibility','off');
grid(ax, 'on');

% X limits and ticks
xlim(ax, [actualDT(max(end-120,1)), fcst_dates(end)+calmonths(1)]);
ticks = datetime(year(actualDT(1)),1,1):calyears(2):fcst_dates(end);
ax.XTick = ticks; ax.XAxis.TickLabelFormat = 'yyyy';

legend(ax,'show','Location','northoutside','Orientation','horizontal'); legend boxoff;
set(ax,'FontSize',13);
 title(ax, sprintf('Inflation (CPI) — last-origin fan (Nov 2025 origin): 1–36 months ahead'));
% last_origin_date = months_origin(end);
% title(ax, sprintf('Inflation (CPI) — last-origin fan (%s origin): 1–36 months ahead', ...
%                   upper(datestr(last_origin_date,'mmm yyyy'))));
hold(ax,'off');

% Save the figure
ensure = @(p) (exist(p,'dir') || mkdir(p));
ensure(fullfile(dropboxFolder,'predictive_densities'));
print(gcf, fullfile(dropboxFolder,'predictive_densities', ...
    'inflation_fanchart_OOS_lastOrigin_36mo.png'), '-dpng', '-r200');



%% === Save last-origin 36-month skew-t parameters (monthly) ===
spec_to_use = 1;

load(fullfile(outputFolder,'sktparam','skewtparam_inflation_OOS.mat'), 'lc_skt','sc_skt','sh_skt','df_skt');
load(fullfile(outputFolder,'actual_inflation_mom_OOS.mat'),'StartEst'); 

% Rebuild months_origin if needed
if ~exist('months_origin','var')
    load(fullfile(outputFolder,'predicted_quantiles_inflation_OOS.mat'),'predicted_quantiles_OOS');
    T_plot = size(predicted_quantiles_OOS,1);
    StartEstDT = datetime(StartEst,'ConvertFrom','datenum');
    months_origin = StartEstDT + calmonths(0:T_plot-1);
end

% Validate horizon availability
H_avail = size(lc_skt,2);
if H_avail < 37
    error('lc_skt has only %d horizons; need >= 37 (0..36). Re-run with horizons = 37.', H_avail);
end

t_last = size(lc_skt, 1);
h_idx36 = 2:37;  % 1..36 ahead

slice3 = @(A) squeeze(A(t_last, h_idx36, min(spec_to_use, size(A,3))));

if ndims(lc_skt) == 3
    lc_last36 = slice3(lc_skt); sc_last36 = slice3(sc_skt);
    sh_last36 = slice3(sh_skt); df_last36 = slice3(df_skt);
else
    lc_last36 = squeeze(lc_skt(t_last, h_idx36));
    sc_last36 = squeeze(sc_skt(t_last, h_idx36));
    sh_last36 = squeeze(sh_skt(t_last, h_idx36));
    df_last36 = squeeze(df_skt(t_last, h_idx36));
end

% As column vectors
lc_last36 = lc_last36(:); sc_last36 = sc_last36(:);
sh_last36 = sh_last36(:); df_last36 = df_last36(:);

% Target dates
last_origin_date = months_origin(t_last);
fcst_dates_36 = (last_origin_date + calmonths(1:36)).';

% Optional: moments
moments_path = fullfile(outputFolder,'sktparam','moments_inflation_OOS.mat');
if exist(moments_path,'file')
    load(moments_path,'fcstmean','fcststdev','fcstskew');
    if ndims(fcstmean) == 3
        fcstmean_last36  = squeeze(fcstmean(t_last,  h_idx36, spec_to_use)).';
        fcststdev_last36 = squeeze(fcststdev(t_last, h_idx36, spec_to_use)).';
        fcstskew_last36  = squeeze(fcstskew(t_last,  h_idx36, spec_to_use)).';
    else
        fcstmean_last36  = squeeze(fcstmean(t_last,  h_idx36)).';
        fcststdev_last36 = squeeze(fcststdev(t_last, h_idx36)).';
        fcstskew_last36  = squeeze(fcstskew(t_last,  h_idx36)).';
    end
else
    fcstmean_last36 = []; fcststdev_last36 = []; fcstskew_last36 = [];
end

% Save monthly (36x1 vectors)
ensure = @(p) (exist(p,'dir') || mkdir(p));
ensure(fullfile(outputFolder,'sktparam'));
out_monthly = fullfile(outputFolder,'sktparam', sprintf('last_origin_skewt_params_36mo_spec%d.mat', spec_to_use));
save(out_monthly, 'lc_last36','sc_last36','sh_last36','df_last36', ...
     'fcstmean_last36','fcststdev_last36','fcstskew_last36', 'fcst_dates_36', ...
     'last_origin_date','spec_to_use');
fprintf('Saved monthly (36) skew-t params to:\n  %s\n', out_monthly);


%% === Build "quarterly" series by taking every 3rd month (3,6,...,36) ===
idx_qtr = 3:3:36;    % 12 values

lc_qtr  = lc_last36(idx_qtr);
sc_qtr  = sc_last36(idx_qtr);
sh_qtr  = sh_last36(idx_qtr);
df_qtr  = df_last36(idx_qtr);

fcst_dates_qtr = fcst_dates_36(idx_qtr);

% (Optional) moments at quarterly cadence
if ~isempty(fcstmean_last36)
    fcstmean_qtr  = fcstmean_last36(idx_qtr);
    fcststdev_qtr = fcststdev_last36(idx_qtr);
    fcstskew_qtr  = fcstskew_last36(idx_qtr);
else
    fcstmean_qtr = []; fcststdev_qtr = []; fcstskew_qtr = [];
end

out_qtr = fullfile(outputFolder,'sktparam', sprintf('last_origin_skewt_params_3y_quarterly_spec%d.mat', spec_to_use));
save(out_qtr, 'lc_qtr','sc_qtr','sh_qtr','df_qtr', ...
     'fcstmean_qtr','fcststdev_qtr','fcstskew_qtr', 'fcst_dates_qtr', ...
     'last_origin_date','spec_to_use','idx_qtr');
fprintf('Saved "quarterly" (every 3rd month from origin) skew-t params to:\n  %s\n', out_qtr);


% 
% % Calendar quarter ends among the 36 forecast months
% is_q_end = ismember(month(fcst_dates_36), [3,6,9,12]);
% lc_q_cal  = lc_last36(is_q_end);
% sc_q_cal  = sc_last36(is_q_end);
% sh_q_cal  = sh_last36(is_q_end);
% df_q_cal  = df_last36(is_q_end);
% fcst_dates_q_cal = fcst_dates_36(is_q_end);
% 
% if ~isempty(fcstmean_last36)
%     fcstmean_q_cal  = fcstmean_last36(is_q_end);
%     fcststdev_q_cal = fcststdev_last36(is_q_end);
%     fcstskew_q_cal  = fcstskew_last36(is_q_end);
% else
%     fcstmean_q_cal = []; fcststdev_q_cal = []; fcstskew_q_cal = [];
% end
% 
% out_q_cal = fullfile(outputFolder,'sktparam', sprintf('last_origin_skewt_params_3y_calendar_quarters_spec%d.mat', spec_to_use));
% save(out_q_cal, 'lc_q_cal','sc_q_cal','sh_q_cal','df_q_cal', ...
%      'fcstmean_q_cal','fcststdev_q_cal','fcstskew_q_cal', 'fcst_dates_q_cal', ...
%      'last_origin_date','spec_to_use');
% fprintf('Saved calendar-aligned quarter-end skew-t params to:\n  %s\n', out_q_cal);
% 


%% === Save quarterly skew-t parameters ONLY (FOR GIULIA) ===

spec_to_use = 1;

% Path to the previously generated quarterly data
inQuarter = fullfile(outputFolder,'sktparam', ...
    sprintf('last_origin_skewt_params_3y_quarterly_spec%d.mat', spec_to_use));

% Load it
S = load(inQuarter);

% Extract ONLY what GIULIA wants
lc_qtr         = S.lc_qtr;
sc_qtr         = S.sc_qtr;
sh_qtr         = S.sh_qtr;
df_qtr         = S.df_qtr;
fcst_dates_qtr = S.fcst_dates_qtr;

% Build a neat output file
outFile = fullfile(outputFolder,'sktparam', ...
    sprintf('quarterly_skewt_params_forSharing_spec%d.mat', spec_to_use));

save(outFile, ...
     'lc_qtr','sc_qtr','sh_qtr','df_qtr','fcst_dates_qtr', ...
     '-v7.3');   % v7.3 ensures compatibility with Python, R, OCTAVE, etc.

fprintf('\nSaved clean 1x12 quarterly .mat file for sharing:\n   %s\n', outFile);

