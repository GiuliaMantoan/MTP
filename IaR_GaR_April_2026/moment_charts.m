%% moment_charts.m
% Charts of forecast moments (mean, standard deviation, skewness) for
% growth and inflation, comparing the BoE fan chart, the semiparametric
% model and the skew-t model.
%
% Unlike the earlier version of this script, growth and inflation are NOT
% assumed to share a horizon frequency:
%   - BoE (growth & inflation) and the growth models (skew-t, semi-param)
%     report horizons h_0..h_12 stepping in QUARTERS (0 to 12 quarters
%     ahead).
%   - The inflation models (skew-t, semi-param) report horizons h_0..h_36
%     stepping in MONTHS (0 to 36 months ahead).
% The script converts everything to a common set of target horizons of
% one quarter, one year and two years ahead by dividing the target
% horizon (in months) by each source's own step size in months, so the
% comparison is apples-to-apples regardless of the underlying frequency.
%
% Series across the three models are aligned by their forecast-origin
% date (inner join), since the BoE, skew-t and semi-param files no longer
% necessarily share the same date range or number of observations.
%
% Covid treatment: the x-axis here is the forecast-ORIGIN date (unlike
% the rolling fan charts, which plot against the target date), so the
% exclusion window does not need to be shifted per horizon. Any series
% observation whose origin falls inside [shadeStart, shadeEnd] is set to
% NaN before plotting, so the model lines show a genuine gap over the
% Covid period (and its base-effect/lag buffer) instead of interpolating
% straight through it. The grey patch is retained purely as a visual
% marker of that same window.

close all; clear; clc;

%% ---- Paths -------------------------------------------------------
basePath   = fileparts(mfilename('fullpath'));   % folder this script lives in
outputsDir = fullfile(basePath, 'Outputs');
outputFolder = fullfile(basePath, 'Moment Charts');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% ---- Style (kept consistent with the previous moment-chart script) ----
models      = {'Fan Charts', 'Semiparam', 'Skew-t'};
colors      = {'b', 'r', 'k'};              % blue, red, black
linestyles  = {'-', '--', ':'};

% Covid exclusion window (origin-date basis). shadeEnd already runs well
% past the acute 2020 shock to absorb the lagged/base-effect distortion
% in the data, mirroring covidEnd + covidLagBuffer in the fan chart script.
shadeStart  = datetime(2020,4,1);
shadeEnd    = datetime(2022,6,30);

metricSheets = {'fcstmean', 'fcststdev', 'fcstskew'};
metricLabels = {'Mean', 'Std. Dev.', 'Skew'};
metricTags   = {'mean', 'std', 'skew'};

targetMonths  = [3, 12, 24, 36];                          % 1Q, 1Y, 2Y, 3Y ahead
horizonTags   = {'1q', '1y', '2y', '3y'};
horizonTitles = {'One Quarter Ahead', 'One Year Ahead', 'Two Years Ahead', 'Three Years Ahead'};

%% ---- Variable configuration ---------------------------------------
cfg = struct([]);

cfg(1).name         = 'Growth';
cfg(1).fileTag       = 'growth';
cfg(1).boeFile       = fullfile(basePath, 'boe_growth_moments.xlsx');
cfg(1).smpFile       = fullfile(outputsDir, 'semi-param_best_spec_gdp_OOS_extended.xlsx');
cfg(1).sktFile       = fullfile(outputsDir, 'skewt_best_spec_gdp_OOS_extended.xlsx');
cfg(1).boeStepMonths = 3;   % BoE growth horizons step in quarters
cfg(1).smpStepMonths = 3;   % growth semi-param horizons step in quarters
cfg(1).sktStepMonths = 3;   % growth skew-t horizons step in quarters
% Growth series can occasionally blow up to implausible extremes (e.g. a
% wide skew-t tail during a volatile quarter), which stretches the
% y-axis so much that the rest of the series looks flat. Clamp the axis
% for growth only, one row per metric (mean/std/skew), matching the
% bounds used previously for the growth charts.
cfg(1).yClamp = [-10 10; -8 8; -4 4];
cfg(1).covidBlank   = true;   % blank moments during Covid, as in the fan charts

cfg(2).name         = 'Inflation';
cfg(2).fileTag       = 'inflation';
cfg(2).boeFile       = fullfile(basePath, 'boe_inflation_moments.xlsx');
cfg(2).smpFile       = fullfile(outputsDir, 'semi-param_best_spec_inflation.xlsx');
cfg(2).sktFile       = fullfile(outputsDir, 'skewt_best_spec_inflation.xlsx');
cfg(2).boeStepMonths = 3;   % BoE inflation horizons step in quarters
cfg(2).smpStepMonths = 1;   % inflation semi-param horizons step in MONTHS
cfg(2).sktStepMonths = 1;   % inflation skew-t horizons step in MONTHS
cfg(2).yClamp = [];         % no axis clamping for inflation
cfg(2).covidBlank   = false;  % leave inflation series intact through Covid (shading only)

%% ---- Main loop over variables --------------------------------------
for v = 1:numel(cfg)

    vcfg = cfg(v);
    fprintf('Processing %s...\n', vcfg.name);

    for m = 1:numel(metricSheets)

        sheetName = metricSheets{m};

        boeTbl = readtable(vcfg.boeFile, 'Sheet', sheetName);
        smpTbl = readtable(vcfg.smpFile, 'Sheet', sheetName);
        sktTbl = readtable(vcfg.sktFile, 'Sheet', sheetName);

        boeDates = boeTbl{:,1};
        smpDates = smpTbl{:,1};
        sktDates = sktTbl{:,1};
        if ~isdatetime(boeDates), boeDates = datetime(boeDates, 'ConvertFrom', 'excel'); end
        if ~isdatetime(smpDates), smpDates = datetime(smpDates, 'ConvertFrom', 'excel'); end
        if ~isdatetime(sktDates), sktDates = datetime(sktDates, 'ConvertFrom', 'excel'); end

        % Align the three sources on their common forecast-origin
        % month. Exact day-of-month differs between sources (e.g. BoE
        % uses quarter-end dates, the model outputs use the 1st of the
        % month), so match on (year, month) rather than the exact date.
        boeKey = dateshift(boeDates, 'start', 'month');
        smpKey = dateshift(smpDates, 'start', 'month');
        sktKey = dateshift(sktDates, 'start', 'month');

        commonDates = intersect(intersect(boeKey, smpKey), sktKey);
        commonDates = sort(commonDates);

        [~, boeRows] = ismember(commonDates, boeKey);
        [~, smpRows] = ismember(commonDates, smpKey);
        [~, sktRows] = ismember(commonDates, sktKey);

        % Covid mask on the (shared) origin-date axis. Computed once per
        % variable/metric since it does not depend on horizon here. Only
        % applied for variables with covidBlank = true (growth); for
        % inflation this stays all-false, so nothing gets blanked and the
        % series is plotted exactly as before, with the grey patch as
        % the only Covid marker.
        if vcfg.covidBlank
            covidMask = commonDates >= shadeStart & commonDates <= shadeEnd;
        else
            covidMask = false(size(commonDates));
        end

        for h = 1:numel(targetMonths)

            tMonths = targetMonths(h);

            boeCol = colIndexForHorizon(tMonths, vcfg.boeStepMonths);
            smpCol = colIndexForHorizon(tMonths, vcfg.smpStepMonths);
            sktCol = colIndexForHorizon(tMonths, vcfg.sktStepMonths);

            boeSeries = boeTbl{boeRows, boeCol};
            smpSeries = smpTbl{smpRows, smpCol};
            sktSeries = sktTbl{sktRows, sktCol};

            % Blank out the Covid period: NaNs naturally break each line
            % over the excluded window instead of interpolating through
            % it, so the chart doesn't imply a stable forecast density
            % existed during that stretch.
            boeSeries(covidMask) = NaN;
            smpSeries(covidMask) = NaN;
            sktSeries(covidMask) = NaN;

            figure('Visible', 'off');
            hold on;
            plot(commonDates, boeSeries, linestyles{1}, 'LineWidth', 1.5, 'Color', colors{1});
            plot(commonDates, smpSeries, linestyles{2}, 'LineWidth', 1.5, 'Color', colors{2});
            plot(commonDates, sktSeries, linestyles{3}, 'LineWidth', 1.5, 'Color', colors{3});
            hold off;

            legend(models, 'Location', 'northwest');
            xlabel('Date');
            ylabel(sprintf('%s (%s)', metricLabels{m}, vcfg.name));
            grid on;

            % Clamp the y-axis for growth (only) so an occasional extreme
            % value doesn't wash out the rest of the series.
            if ~isempty(vcfg.yClamp)
                yl = ylim;
                yl(1) = max(yl(1), vcfg.yClamp(m,1));
                yl(2) = min(yl(2), vcfg.yClamp(m,2));
                ylim(yl);
            end

            yl = ylim;
            hPatch = patch( ...
                [shadeStart shadeStart shadeEnd shadeEnd], ...
                [yl(1) yl(2) yl(2) yl(1)], ...
                [0.6 0.6 0.6], ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.3);
            uistack(hPatch, 'bottom');
            hPatch.HandleVisibility = 'off';
            grid on;
            yline(0, 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off');

            % title(sprintf('%s: %s', vcfg.name, horizonTitles{h}));

            fname = sprintf('%s_%s_%s.png', vcfg.fileTag, metricTags{m}, horizonTags{h});
            saveas(gcf, fullfile(outputFolder, fname));
            close(gcf);
        end
    end
end

fprintf('Done. Charts saved to: %s\n', outputFolder);

%% ---- Local functions -------------------------------------------------
function colIdx = colIndexForHorizon(targetMonths, stepMonths)
% Maps a target horizon (in months) to the column index within a moments
% table (column 1 is Dates, column 2 is h_0, column 3 is h_1, ...), given
% that particular source's horizon step size (in months).
    nSteps = targetMonths / stepMonths;
    if mod(nSteps, 1) ~= 0
        error('Target horizon of %d months is not reachable with a step size of %d months.', ...
            targetMonths, stepMonths);
    end
    colIdx = nSteps + 2; % +1 for h_0 offset, +1 for the Dates column
end