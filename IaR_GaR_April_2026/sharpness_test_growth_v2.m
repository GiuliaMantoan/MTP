%% Sharpness test %%

% This code executes a width and graphic based sharpness evaluation
% Authors: David Aikman, Rhys Bidder, Simon Lloyd, Giulia Mantoan, Simone
% Maso, Aditya Mori and Matthew Tong

%% Loading the pdf values.

% This code does not require loading the actual data or the creation of
% PITs but simply the densities.

% This code focuses on the Bank of England Fan Charts and the Skew-t
% density forecasts.

close all; clear; clc;

%% Settings

set(0,'defaultAxesFontName', 'Times'); % font for chart axis
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1) % line style
rng('default'); % defaul random number generateor
rng(0);

cd 'C:\Users\adity\Dropbox\BoE-KCL Macro Forecasting\' % chage here the cd
addpath('Data\') % data input
addpath('Codes\intermediate_codes\')  % intemediate code (data loading..)
addpath('Codes\functions\')
addpath('Codes\functions\azzalini')  % skew t functions

% Dates for importing relevant data
start_date = datetime('30-Sep-2004', 'InputFormat', 'dd-MMM-yyyy'); % first forecasted period (Q3-2004)
end_date   = datetime('31-Mar-2024', 'InputFormat', 'dd-MMM-yyyy'); % last date with actual data
end_fcst   = datetime('30-Mar-2027', 'InputFormat', 'dd-MMM-yyyy'); % forecast end date
covid_date = datetime('30-Jun-2020', 'InputFormat', 'dd-MMM-yyyy');
covid_end_date = datetime('30-Jun-2022', 'InputFormat', 'dd-MMM-yyyy');

%% Constructing the Fan Charts

% Load the BoE forecast data
boefsctdata_growth;

outDir = 'C:\Users\adity\Dropbox\BoE-KCL Macro Forecasting\Codes\Auxiliary and trial codes AM\Sharpness Charts';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% At each forecast horizon and forecast date combination, need to evaluate the distribution
% Generating the support

nSupport = 500;

% compute which quarter index is COVID (Q2-2020)
covid = (year(covid_date) - year(start_date)) * 4 + (ceil(month(covid_date)/3) - ceil(month(start_date)/3));
covid_end = (year(covid_end_date) - year(start_date)) * 4 + (ceil(month(covid_end_date)/3) - ceil(month(start_date)/3));

modevar = mtestdata(:,1:3:end);
modevar(:, covid+1) = 1.5; % fake data because none exists for this date
meanvar = mtestdata(:,2:3:end);
meanvar(:, covid+1) = 1.4; % fake data
dispvar = mtestdata(:,3:3:end);
dispvar(:, covid+1) = 0.3; % fake data

[H, T] = size(meanvar);

fan_pdf      = zeros(H, T, nSupport);
support_grid = zeros(H, T, nSupport);

for i = 1:H
    for j = 1:T
        if meanvar(i,j) - modevar(i,j) ~= 0
            gam = stdtogam(meanvar(i,j), modevar(i,j), dispvar(i,j)); % Skewness in  Britton, Fisher and Whitley formulation
            sig = mom2g(dispvar(i,j), gam); % get the variance in the common formulation
            [m, s1, s2] = momtopar(meanvar(i,j), modevar(i,j), sig); % get lhs and rhs std dev (m is just the mode)
            xlow = m - 4*s1;
            xhigh = m + 4*s2;
            y = linspace(xlow, xhigh, nSupport);
            fan_pdf(i, j, :) = ftp(y, m, s1, s2);
            support_grid(i,j,:)   = y;
        else
            xlow = meanvar(i,j) - 4*dispvar(i,j);
            xhigh = meanvar(i,j) + 4*dispvar(i,j);
            y = linspace(xlow, xhigh, nSupport);
            fan_pdf(i, j, :) = normpdf(y, meanvar(i, j), dispvar(i,j)); % standardized value for growth if mean = mode
            support_grid(i,j,:)   = y;
        end
    end
end

% Now need to remove the pdf that is affected by covid

for tIdx = covid+1 : covid_end+1
    startCol = tIdx - 12;        % the first column from which we begin deleting
    for k = 0 : 12
        col = startCol + k;      % current column
        row = H - k;             % which row to delete (k=0→last row, k=12→row  T_t-12)
        % sanity check
        if col>=1 && col<=T  &&  row>=1 && row<=H
            fan_pdf(row, col, :) = NaN;
        end
    end
end

fan_pdf(:, (covid+1):(covid+5), :) = NaN; % previously was fan_pdf(:, covid+1, :) = NaN; this has been changed to deal with the fact that the dispersion variable causes the fan to be particularly wide. 

clear tIdx covid covid_end k col row

%% Skew t pdf construction

load("C:\Users\adity\Dropbox\BoE-KCL Macro Forecasting\Outputs\sktparam\skewtparam_all_spec_gdp_yoy_4_cat.mat")
start_date = datetime('31-Mar-2004', 'InputFormat', 'dd-MMM-yyyy');
end_date   = datetime('31-Mar-2025', 'InputFormat', 'dd-MMM-yyyy'); % last date with actual data.
end_fcst   = datetime('31-Mar-2028', 'InputFormat', 'dd-MMM-yyyy'); % forecast end date

% Covid treatment
% the covid and covid_end index will be slightly different given different starting
% dates
covid = (year(covid_date) - year(start_date)) * 4 + (ceil(month(covid_date)/3) - ceil(month(start_date)/3));
covid_end = (year(covid_end_date) - year(start_date)) * 4 + (ceil(month(covid_end_date)/3) - ceil(month(start_date)/3));

% List of variable base names
varnames = {'lc_skt', 'sc_skt', 'sh_skt', 'df_skt'};

% Loop through and remove the covid+1 column if it exists
for i = 1:length(varnames)
    varname = varnames{i};
    var = eval(varname);
    var = squeeze(var(:, :, 2)); % takes the second specification from the four variables specification.
    assignin('base', varname, var);
end

skewt_pdf = zeros(T, H, nSupport);
support_grid_skewt = zeros(T, H, nSupport);

for i = 1:T
    for j = 1:H
        df = df_skt(i,j);
        location = lc_skt(i,j);
        scale = sc_skt(i,j);
        shape = sh_skt(i,j);
        xlow  = lc_skt(i, j) - 4*sc_skt(i, j);
        xhigh  = lc_skt(i, j) + 4*sc_skt(i, j);
        x = linspace(xlow, xhigh, nSupport);
        skewt_pdf(i,j, :) = dskt(x,location,scale,shape,df);
        support_grid_skewt(i, j, :) = x;
    end
end

for tIdx = covid+1 : covid_end+1
    startCol = tIdx - 12;        % the first column from which we begin deleting
    for k = 0 : 12
        col = startCol + k;      % current column
        row = H - k;             % which row to delete (k=0→last row, k=12→row  T_t-12)
        % sanity check
        if col>=1 && col<=T  &&  row>=1 && row<=H
            skewt_pdf(col, row, :) = NaN;
        end
    end
end

skewt_pdf((covid+1):(covid+5), :, :) = NaN;

% Finding the widths

% fans width initialising
width50_fan  = zeros(H, T);
width90_fan  = zeros(H, T);
widthb20_fan = zeros(H, T);
widtht20_fan = zeros(H, T);
width30_fan  = zeros(H, T);

% Initialising fan quantiles
q05_fan = zeros(H, T);
q10_fan = zeros(H, T);
q25_fan = zeros(H, T);
q50_fan = zeros(H, T);
q75_fan = zeros(H, T);
q90_fan = zeros(H, T);
q95_fan = zeros(H, T);

% skew-t width initialising
width50_skt  = zeros(T, H);
width90_skt  = zeros(T, H);
widthb20_skt = zeros(T, H);
widtht20_skt = zeros(T, H);
width30_skt  = zeros(T, H);

% Initialising skew-t quantiles
q05_skt = zeros(T, H);
q10_skt = zeros(T, H);
q25_skt = zeros(T, H);
q50_skt = zeros(T, H);
q75_skt = zeros(T, H);
q90_skt = zeros(T, H);
q95_skt = zeros(T, H);


%% Fan-chart widths
dx_fan = support_grid(1,1,2) - support_grid(1,1,1);

for i = 1:H
    for j = 1:T
        x_ij = squeeze(support_grid(i,j,:));
        p_ij = squeeze(fan_pdf(i,j,:));

        [cu, xu] = build_strict_cdf(x_ij, p_ij);
        if isempty(cu)
            width50_fan(i,j)=NaN; width30_fan(i,j)=NaN; width90_fan(i,j)=NaN;
            widthb20_fan(i,j)=NaN; widtht20_fan(i,j)=NaN;
            q05_fan(i,j)=NaN; q10_fan(i,j)=NaN; q25_fan(i,j)=NaN; q50_fan(i,j)=NaN;
            q75_fan(i,j)=NaN; q90_fan(i,j)=NaN; q95_fan(i,j)=NaN;
            continue
        end

        q05 = safe_quantile(cu, xu, 0.05);
        q10 = safe_quantile(cu, xu, 0.10);
        q25 = safe_quantile(cu, xu, 0.25);
        q35 = safe_quantile(cu, xu, 0.35);
        q50 = safe_quantile(cu, xu, 0.50);
        q65 = safe_quantile(cu, xu, 0.65);
        q75 = safe_quantile(cu, xu, 0.75);
        q90 = safe_quantile(cu, xu, 0.90);
        q95 = safe_quantile(cu, xu, 0.95);

        width50_fan(i,j)  = q75 - q25;
        width30_fan(i,j)  = q65 - q35;
        width90_fan(i,j)  = q95 - q05;
        widthb20_fan(i,j) = q25 - q05;
        widtht20_fan(i,j) = q95 - q75;

        q05_fan(i,j)=q05; q10_fan(i,j)=q10; q25_fan(i,j)=q25; q50_fan(i,j)=q50;
        q75_fan(i,j)=q75; q90_fan(i,j)=q90; q95_fan(i,j)=q95;

    end
end

% fan chart quantile means for plotting across forecast dates
q05_fan_mean = mean(q05_fan, 2, 'omitnan');
q10_fan_mean = mean(q10_fan, 2, 'omitnan');
q25_fan_mean = mean(q25_fan, 2, 'omitnan');
q50_fan_mean = mean(q50_fan, 2, 'omitnan');
q75_fan_mean = mean(q75_fan, 2, 'omitnan');
q90_fan_mean = mean(q90_fan, 2, 'omitnan');
q95_fan_mean = mean(q95_fan, 2, 'omitnan');

clear q05 q10 q25 q35 q65 q75 q90 q95 x_ij cu xu p_ij

%% Skew-t widths
dx_skt = support_grid_skewt(1,1,2) - support_grid_skewt(1,1,1);

for i = 1:T
    for j = 1:H
        x_ij = squeeze(support_grid_skewt(i,j,:));
        p_ij = squeeze(skewt_pdf(i,j,:));

        [cu, xu] = build_strict_cdf(x_ij, p_ij);
        if isempty(cu)
            width50_skt(i,j)=NaN; width30_skt(i,j)=NaN; width90_skt(i,j)=NaN;
            widthb20_skt(i,j)=NaN; widtht20_skt(i,j)=NaN;
            q05_skt(i,j)=NaN; q10_skt(i,j)=NaN; q25_skt(i,j)=NaN; q50_skt(i,j)=NaN;
            q75_skt(i,j)=NaN; q90_skt(i,j)=NaN; q95_skt(i,j)=NaN;
            continue
        end

        q05 = safe_quantile(cu, xu, 0.05);
        q10 = safe_quantile(cu, xu, 0.10);
        q25 = safe_quantile(cu, xu, 0.25);
        q35 = safe_quantile(cu, xu, 0.35);
        q50 = safe_quantile(cu, xu, 0.50);
        q65 = safe_quantile(cu, xu, 0.65);
        q75 = safe_quantile(cu, xu, 0.75);
        q90 = safe_quantile(cu, xu, 0.90);
        q95 = safe_quantile(cu, xu, 0.95);

        width50_skt(i,j)  = q75 - q25;
        width30_skt(i,j)  = q65 - q35;
        width90_skt(i,j)  = q95 - q05;
        widthb20_skt(i,j) = q25 - q05;
        widtht20_skt(i,j) = q95 - q75;

        q05_skt(i,j)=q05; q10_skt(i,j)=q10; q25_skt(i,j)=q25; q50_skt(i,j)=q50;
        q75_skt(i,j)=q75; q90_skt(i,j)=q90; q95_skt(i,j)=q95;

    end
end

% Skew-t quantiles means across forecast dates
q05_skt_mean = mean(q05_skt, 1, 'omitnan');
q10_skt_mean = mean(q10_skt, 1, 'omitnan');
q25_skt_mean = mean(q25_skt, 1, 'omitnan');
q50_skt_mean = mean(q50_skt, 1, 'omitnan');
q75_skt_mean = mean(q75_skt, 1, 'omitnan');
q90_skt_mean = mean(q90_skt, 1, 'omitnan');
q95_skt_mean = mean(q95_skt, 1, 'omitnan');

% Mean width of the fan charts around specific intervals
fan_width_mean_50 = mean(width50_fan, 2, 'omitnan');
fan_width_mean_90 = mean(width90_fan, 2, 'omitnan');
fan_width_mean_30 = mean(width30_fan, 2, 'omitnan');
fan_width_mean_b20= mean(widthb20_fan, 2, 'omitnan');
fan_width_mean_t20= mean(widtht20_fan, 2, 'omitnan');

% Mean width of the skew-t densities around specific intervals
skewt_width_mean_50 = mean(width50_skt, 1, 'omitnan');
skewt_width_mean_90 = mean(width90_skt, 1, 'omitnan');
skewt_width_mean_30 = mean(width30_skt, 1, 'omitnan');
skewt_width_mean_b20= mean(widthb20_skt, 1, 'omitnan');
skewt_width_mean_t20= mean(widtht20_skt, 1, 'omitnan');

fan_mu50 = median(width50_fan, 2, 'omitnan');
fan_sd50 = std(width50_fan, 0, 2, 'omitnan');
fan_N50  = sum(isfinite(width50_fan), 2);
fan_se50 = fan_sd50 ./ sqrt(fan_N50);   % if you prefer s.e. bands

fan_mu30 = median(width30_fan, 2, 'omitnan');
fan_sd30 = std(width30_fan, 0, 2, 'omitnan'); fan_N30 = sum(isfinite(width30_fan), 2); fan_se30 = fan_sd30./sqrt(fan_N30);

fan_mu90 = median(width90_fan, 2, 'omitnan');
fan_sd90 = std(width90_fan, 0, 2, 'omitnan'); fan_N90 = sum(isfinite(width90_fan), 2); fan_se90 = fan_sd90./sqrt(fan_N90);

fan_mub20 = median(widthb20_fan, 2, 'omitnan');
fan_sdb20 = std(widthb20_fan, 0, 2, 'omitnan'); fan_Nb20 = sum(isfinite(widthb20_fan), 2); fan_seb20 = fan_sdb20./sqrt(fan_Nb20);

fan_mut20 = median(widtht20_fan, 2, 'omitnan');
fan_sdt20 = std(widtht20_fan, 0, 2, 'omitnan'); fan_Nt20 = sum(isfinite(widtht20_fan), 2); fan_set20 = fan_sdt20./sqrt(fan_Nt20);

% Skew-t (T×H → along dim=1). Transpose to column vectors for plotting
skt_mu50 = median(width50_skt, 1, 'omitnan')'; 
skt_sd50 = std(width50_skt, 0, 1, 'omitnan')';
skt_N50  = sum(isfinite(width50_skt), 1)'; 
skt_se50 = skt_sd50 ./ sqrt(skt_N50);

skt_mu30 = median(width30_skt, 1, 'omitnan')'; skt_sd30 = std(width30_skt, 0, 1, 'omitnan')'; skt_N30 = sum(isfinite(width30_skt),1)'; skt_se30 = skt_sd30./sqrt(skt_N30);
skt_mu90 = median(width90_skt, 1, 'omitnan')'; skt_sd90 = std(width90_skt, 0, 1, 'omitnan')'; skt_N90 = sum(isfinite(width90_skt),1)'; skt_se90 = skt_sd90./sqrt(skt_N90);
skt_mub20= median(widthb20_skt,1, 'omitnan')'; skt_sdb20= std(widthb20_skt,0, 1, 'omitnan')'; skt_Nb20= sum(isfinite(widthb20_skt),1)'; skt_seb20= skt_sdb20./sqrt(skt_Nb20);
skt_mut20= median(widtht20_skt,1, 'omitnan')'; skt_sdt20= std(widtht20_skt,0, 1, 'omitnan')'; skt_Nt20= sum(isfinite(widtht20_skt),1)'; skt_set20= skt_sdt20./sqrt(skt_Nt20);

horizons = (1:numel(fan_mu50))';

%% PLOT: CENTRAL 30% 
figure; hold on; box on; grid on;

% Fan: mean ± 1*sd (use fan_se50 instead if you want s.e. bands)
plot_with_band(horizons, fan_mu30, fan_sd30, [0.7 0.8 1.0], 'b-', 'Fan Median Width');

% Skew-t
plot_with_band(horizons, skt_mu30, skt_sd30, [1.0 0.8 0.8], 'r--', 'Skew-t Median Width');

xlabel('Horizon'); ylabel('Width');
title('Median Width in Central 30% Interval (±1 SD)');
legend('Location','southoutside','Orientation','horizontal'); 
exportgraphics(gcf, fullfile(outDir, 'width30_growth_bands.png'), 'Resolution', 300);


%% PLOT: 50% CENTRAL

figure; hold on; box on; grid on;

% Fan: median ± 1*sd (use fan_se50 instead if you want s.e. bands)
plot_with_band(horizons, fan_mu50, fan_sd50, [0.7 0.8 1.0], 'b-', 'Fan Median Width');

% Skew-t
plot_with_band(horizons, skt_mu50, skt_sd50, [1.0 0.8 0.8], 'r--', 'Skew-t Median Width');

xlabel('Horizon'); ylabel('Width');
title('Median Width in Central 50% Interval (±1 SD)');
legend('Location','southoutside','Orientation','horizontal'); 
exportgraphics(gcf, fullfile(outDir, 'width50_growth_bands.png'), 'Resolution', 300);

%% PLOT: 90% CENTRAL

figure; hold on; box on; grid on;

% Fan: mean ± 1*sd (use fan_se50 instead if you want s.e. bands)
plot_with_band(horizons, fan_mu90, fan_sd90, [0.7 0.8 1.0], 'b-', 'Fan Median Width');

% Skew-t
plot_with_band(horizons, skt_mu90, skt_sd90, [1.0 0.8 0.8], 'r--', 'Skew-t Median Width');

xlabel('Horizon'); ylabel('Width');
title('Median Width in Central 90% Interval (±1 SD)');
legend('Location','southoutside','Orientation','horizontal'); 
exportgraphics(gcf, fullfile(outDir, 'width90_growth_bands.png'), 'Resolution', 300);

%% PLOT: BOTTOM 20%

figure; hold on; box on; grid on;

% Fan: mean ± 1*sd (use fan_se50 instead if you want s.e. bands)
plot_with_band(horizons, fan_mub20, fan_sdb20, [0.7 0.8 1.0], 'b-', 'Fan Median Width');

% Skew-t
plot_with_band(horizons, skt_mub20, skt_sdb20, [1.0 0.8 0.8], 'r--', 'Skew-t Median Width');

xlabel('Horizon'); ylabel('Width');
title('Median Width in Bottom 20% Interval (±1 SD)');
legend('Location','southoutside','Orientation','horizontal'); 
exportgraphics(gcf, fullfile(outDir, 'widthb20_growth_bands.png'), 'Resolution', 300);

%% PLOT: TOP 20%

figure; hold on; box on; grid on;

% Fan: mean ± 1*sd (use fan_se50 instead if you want s.e. bands)
plot_with_band(horizons, fan_mut20, fan_sdt20, [0.7 0.8 1.0], 'b-', 'Fan Median Width');

% Skew-t
plot_with_band(horizons, skt_mut20, skt_sdt20, [1.0 0.8 0.8], 'r--', 'Skew-t Median Width');

xlabel('Horizon'); ylabel('Width');
title('Median Width in Top 20% Interval (±1 SD)');
legend('Location','southoutside','Orientation','horizontal'); 
exportgraphics(gcf, fullfile(outDir, 'widtht20_growth_bands.png'), 'Resolution', 300);

%% Composite (2×2) figure: Fan vs Skew-t width comparisons
% Choose a descriptive label for the export/figure title
labelStr = 'Growth Forecast Sharpness Comparison';   % <- change if you want a different label

% Align horizons in case fan/skew-t lengths differ
H_use = min([numel(horizons), numel(skt_mu50), numel(fan_mu50), ...
             numel(fan_mu30), numel(skt_mu30), ...
             numel(fan_mub20), numel(skt_mub20), ...
             numel(fan_mut20), numel(skt_mut20)]);
ix    = 1:H_use;
x     = horizons(ix);

% Pack the four panels
panels = {
    struct('title','Central 30%','fan_mu',fan_mu30(ix),'fan_sd',fan_sd30(ix),'skt_mu',skt_mu30(ix),'skt_sd',skt_sd30(ix))
    struct('title','Central 50%','fan_mu',fan_mu50(ix),'fan_sd',fan_sd50(ix),'skt_mu',skt_mu50(ix),'skt_sd',skt_sd50(ix))
    struct('title','Bottom 20%','fan_mu',fan_mub20(ix),'fan_sd',fan_sdb20(ix),'skt_mu',skt_mub20(ix),'skt_sd',skt_sdb20(ix))
    struct('title','Top  20%','fan_mu',fan_mut20(ix),'fan_sd',fan_sdt20(ix),'skt_mu',skt_mut20(ix),'skt_sd',skt_sdt20(ix))
    };

% Compute common y-lims across panels (nice for visual comparability)
ymin = inf; ymax = -inf;
for p = 1:numel(panels)
    lo = [panels{p}.fan_mu - panels{p}.fan_sd; panels{p}.skt_mu - panels{p}.skt_sd];
    hi = [panels{p}.fan_mu + panels{p}.fan_sd; panels{p}.skt_mu + panels{p}.skt_sd];
    ymin = min([ymin; lo], [], 'omitnan');
    ymax = max([ymax; hi], [], 'omitnan');
end
if ~isfinite(ymin) || ~isfinite(ymax)
    ymin = 0; ymax = 1; % fallback
end

fig = figure('Color','w','Position',[100 100 1200 800]);
t   = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','loose');

% Title (sgtitle or title(t,...) depending on MATLAB version)
try
    sgtitle(t, [labelStr ' : Width Comparison']);
catch
    title(t, [labelStr ' : Width Comparison']);
end

legendHandles = gobjects(1,2); % capture from first panel for a global legend

for p = 1:numel(panels)
    ax = nexttile(t,p); hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    % Fan (blue) band + line
    hFan = plot_with_band_ax(ax, x, panels{p}.fan_mu, panels{p}.fan_sd, [0.7 0.8 1.0], 'b-', 'Fan Median Width');

    % Skew-t (red dashed) band + line
    hSkt = plot_with_band_ax(ax, x, panels{p}.skt_mu, panels{p}.skt_sd, [1.0 0.8 0.8], 'r--', 'Skew-t Median Width');

    xlabel(ax,'Horizon'); ylabel(ax,'Width');
    title(ax, panels{p}.title);
    ylim(ax, [ymin ymax]);

    if p == 1
        legendHandles(1) = hFan;  % capture the LINE handles (not the patches)
        legendHandles(2) = hSkt;
    end
end

% Global legend under the grid
lgd = legend(legendHandles, {'Fan Median Width','Skew-t Median Width'}, ...
    'Orientation','horizontal','NumColumns',2);
try
    lgd.Layout.Tile = 'south';  % requires R2020b+
catch
    % Fallback: place legend manually at bottom center if Layout.Tile unsupported
    set(lgd,'Position',[0.35 0.01 0.3 0.05]);
end

% Save once
exportgraphics(fig, fullfile(outDir, ['width_composite_' regexprep(labelStr,'\s+','_') '.png']), 'Resolution', 300);
% close(fig);  % optional

function hLine = plot_with_band_ax(ax, x, y, yerr, faceColor, lineSpec, labelName)
    % shaded band (hidden from legend)
    xx = [x; flipud(x)];
    yy = [y - yerr; flipud(y + yerr)];
    patch('Parent',ax,'XData',xx,'YData',yy,'FaceColor',faceColor, ...
          'EdgeColor','none','FaceAlpha',0.5,'HandleVisibility','off');
    % main line (shows in legend)
    hLine = plot(ax, x, y, lineSpec, 'LineWidth',1.5, 'DisplayName',labelName);
end

function plot_with_band(x, y, yerr, faceColor, lineSpec, labelName)
    % shaded area
    fill([x; flipud(x)], ...
         [y - yerr; flipud(y + yerr)], ...
         faceColor, ...
         'EdgeColor','none', ...
         'FaceAlpha',0.5, ...
         'HandleVisibility','off'); % <- this hides the patch from legend

    % main line (this is what shows in legend)
    plot(x, y, lineSpec, 'LineWidth',1.5, 'DisplayName',labelName);
end

function [cu, xu] = build_strict_cdf(x, p)
% Remove non-finite & zero-length issues, return cu (probabilities) and xu (values)
x = x(:); p = p(:);
valid = isfinite(x) & isfinite(p);
if nnz(valid) < 2
    cu = []; xu = [];
    return
end
x = x(valid);
p = p(valid);

% local step sizes on the filtered support
dx = [diff(x); x(end)-x(end-1)];
% guard against non-positive or non-finite dx (shouldn't happen with linspace, but safe)
good_dx = isfinite(dx) & (dx > 0);
x  = x(good_dx);
p  = p(good_dx);
dx = dx(good_dx);
if numel(x) < 2
    cu = []; xu = [];
    return
end

cdf = cumsum(p .* dx);
S = cdf(end);
if ~isfinite(S) || S <= 0
    cu = []; xu = [];
    return
end
cdf = cdf / S;

% ensure strictly increasing & finite
[cu, ia] = unique(cdf, 'stable');
xu = x(ia);

% final scrub of non-finites
fin = isfinite(cu) & isfinite(xu);
cu = cu(fin); xu = xu(fin);

if numel(cu) < 2
    cu = []; xu = [];
end
end

function q = safe_quantile(cu, xu, p)
% Return NaN if we can't safely interpolate; clamp p into [min(cu), max(cu)]
if isempty(cu) || isempty(xu)
    q = NaN; return
end
p = min(max(p, cu(1)), cu(end));
% 'linear' is safest; if you prefer 'makima', keep it—but linear avoids overshoot
q = interp1(cu, xu, p, 'linear');
end

function draws = inverse_transform_sampling(x, pdf_vals, N)
[cu, xu] = build_strict_cdf(x, pdf_vals);
if isempty(cu)
    draws = nan(N,1); return
end
U = rand(N,1);
U = min(max(U, cu(1)), cu(end));
draws = interp1(cu, xu, U, 'linear');
end
