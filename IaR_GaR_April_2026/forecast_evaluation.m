                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %%% Forecast evaluation code %%%
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Authors: David Aikman, Rhys Bidder, Simone Maso, Aditya Mori 

close all; clear; clc;

%% Settings
 
set(0,'defaultAxesFontName', 'Times'); % font for chart axis 
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1) % line style 
rng('default'); % defaul random number generateor 
rng(0);     

scriptDir    = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'intermediate_codes'));
addpath(fullfile(scriptDir, 'functions'));
addpath(fullfile(scriptDir, 'functions', 'azzalini'));
outputFolder = fullfile(scriptDir, 'Outputs');
evalDir      = fullfile(outputFolder, 'forecast_evaluation');
if ~exist(evalDir, 'dir'), mkdir(evalDir); end

% file names
fullFileName = fullfile(scriptDir, 'IaRDataRaw_monthly_M.xlsx'); % raw data

%% Import Data

%%%%%%%%%%%%%%%%%%%%%%
% Data loading panel %
%%%%%%%%%%%%%%%%%%%%%%

start_date = datetime('31-Mar-2010', 'InputFormat', 'dd-MMM-yyyy'); % first forecasted period (Q3-2004) 
end_date   = datetime('30-Sep-2022', 'InputFormat', 'dd-MMM-yyyy'); % last date with actual data
end_fcst   = datetime('30-Sep-2025', 'InputFormat', 'dd-MMM-yyyy'); % forecast end date
covid_date = datetime('30-Jun-2020', 'InputFormat', 'dd-MMM-yyyy'); 
varnames = {'g4cpi'}; % actual data variables you want to import 
ctrynames = {'UK'}; % actual data country you want to import 
momentlist = {'fcstmean', 'fcststdev', 'fcstskew'}; % moments for OLS and skewt

boefsctdata; % this code focus only on inflation (potentially we can extend it to GDP, unmpl..)

% mtestdata %

% look at the data in block of 3 columns:
% col 1: mode
% col 2: mean
% col 3: sd (sqrt of scale)
% each "block" correspond to a forecast period 
% rows FORECASTED (mean/mode/sd) inflation x-period ahead with x = 1, .. , 9 
% example: (say that the first period is Q3-2004) row 1 refers to the 1-step ahead forecast 
% example: colomn 1 represent the forecasted mode inflation starting
% from Q3-2004 onward. colomn 2 represent the forecasted mean inflation starting
% from Q3-2004 onward...
% colomn 4 represent the forecasted mode inflation starting
% from Q4-2004 onward.

actualdata_monthly;

% actualvar is 37 × n_monthly after actualdata_monthly.
% Reduce to 13 fan-chart horizons × quarter-end origins only (Mar/Jun/Sep/Dec).
% dateVec is preserved by actualdata_monthly so we can use it for filtering.
qtr_mask_eval = ismember(month(dateVec), [3 6 9 12]);
actualvar     = actualvar(1:13, qtr_mask_eval);   % 13 horizons × nQOrig

% actualvar layout:
% rows 1-13 : ACTUAL inflation 0- to 12-quarters ahead
% cols      : one per quarter-end origin (Mar 2010, Jun 2010, …)
 

%% Import the PITs 

%%%%%%%%%%%%%%%%%%%
% Model selection %
%%%%%%%%%%%%%%%%%%%

modelfcst = 0; % 0 is boe, 1 is ols, 3 is skewt
modellist = {'boe','ols', 'qreg', 'skewt'};


if modelfcst == 0

    %%%%%%%%%%%
    %   BOE   % (here we also compute mean, std dev and skew)
    %%%%%%%%%%%

    % Extract the corresponding columns from mominf (mode | mean | dispersion)
    modevar = mtestdata(:, 1:3:end);  % mode 
    meanvar = mtestdata(:, 2:3:end); % mean
    square_root_dispersion_var = mtestdata(:, 3:3:end); % dispersion 
    
    % generate a zero matrix to store the values 
    zinf = NaN(size(actualvar));
    stddevvar = NaN(size(actualvar));
    skewvar = NaN(size(actualvar));

    % loop
    
    for i = 1:size(meanvar, 1)
    
        for j = 1:size(meanvar, 2)
    
            %if ~isnan(actualvar(i, j)) 
    
                if meanvar(i,j) - modevar(i,j) ~= 0
    
                    gam = stdtogam(meanvar(i,j), modevar(i,j), square_root_dispersion_var(i,j)); % Skewness in  Britton, Fisher and Whitley formulation 
                    sig = mom2g(square_root_dispersion_var(i,j), gam); % get the variance in the common formulation
                   
                    % parameter of the distribution 
                    [m, s1, s2] = momtopar(meanvar(i,j), modevar(i,j), sig); % get lhs and rhs std dev (m is just the mode) 
                    
                    if ~isnan(actualvar(i, j)) 
                        zinf(i, j) = integral(@(y) ftp(y, m, s1, s2), -3, actualvar(i, j)); % get the cdf
                    else
                        zinf(i, j) = NaN;
                    end

                    % get sd and skew
                    stddevvar(i,j) = sqrt(sig);
                    term1_skw = sqrt(2/pi) * (s2 - s1);
                    term2_skw = (4/pi - 1) .* ((s2 - s1).^2 + s1 .* s2);
                    skewvar(i,j) =  term1_skw .* term2_skw;

                else
                   
                    if ~isnan(actualvar(i, j)) 
                        zinf(i, j) = normcdf((actualvar(i, j) - meanvar(i,j)) / square_root_dispersion_var(i,j)); % standardized value for inflation if mean = mode 
                    else
                        zinf(i, j) = NaN;
                    end

                    % get sd and skew
                    stddevvar(i,j) = square_root_dispersion_var(i,j);
                    skewvar(i,j) =  0;

                end
           % end
        end
    end

    % make sure no NaN col 
    zinf = zinf(:, ~all(isnan(zinf), 1));

    % save the outcome for Boe Forecast (no covid zero-column for inflation)
    fcstmean  = meanvar;
    fcststdev = stddevvar;
    fcstskew  = skewvar;

    % settings
    nQOrig_boe   = size(fcstmean, 2);
    quarterDates = start_date + calmonths(3*(0:nQOrig_boe-1));
    momentlist = {'fcstmean', 'fcststdev', 'fcstskew'};
    horizons = 13;
    filename = fullfile('boe.xlsx');

    % Loop over each variable in the list
    forecastNames = cellstr(strcat('h_', string(0:horizons-1)));
    for i = 1:length(momentlist)

        % Get the current forecast data  (13 × nCols)
        currentVar = eval(momentlist{i});
        nCols  = size(currentVar, 2);
        qDates = start_date + calmonths(3*(0:nCols-1));   % 1×nCols datetime

        % Build table: nCols rows × (1 date + 13 horizon) cols
        T = array2table(currentVar', 'VariableNames', forecastNames);
        T.Dates = qDates(:);
        T = T(:, [{'Dates'}, forecastNames]);

        writetable(T, filename, 'Sheet', momentlist{i});

    end

    clearvars stddevvar  term1_skw  term2_skw  skewvar  fcstmean fcststdev fcstskew quarterDates  momentlist horizons filename currentVar T forecastNames 

elseif modelfcst == 1
 
    %%%%%%%%%%%
    %   OLS   %
    %%%%%%%%%%%

    % file name 
    filename = 'ols.xlsx';
    
    
    % Loop through momentlist and read each sheet into an array
    for i = 1:length(momentlist)
        sheetname = momentlist{i};
        
        % Read the table from the corresponding sheet
        T = readtable(filename, 'Sheet', sheetname);
        
        % Convert to array (excluding first column, assumed to be dates)
        data_array = table2array(T(:, 2:end))';
        
        % delete covid 
        % Drop columns (63*3 + 1) to (63*3 + 4) of mtestdata if mtestdata has at least 190 columns
        if size(data_array, 2) >= (covid + 1)
            data_array(:, covid + 1) = [];
        end
        
        % Assign dynamically to variable with _array suffix
        assignin('base', [sheetname '_array'], data_array);
    
    end
    
    % generate a zero matrix to store the values 
    zinf = NaN(size(actualvar));
    
    % loop
    
    for i = 1:size(fcstmean_array, 1)
        
        for j = 1:size(fcstmean_array, 2)
        
                if ~isnan(actualvar(i, j)) 
    
                    zinf(i, j) = normcdf(actualvar(i,j),fcstmean_array(i,j),fcststdev_array(i,j)); 
    
                end
        end
    end

elseif modelfcst == 3 

    %%%%%%%%%%%
    %  SKEWT  %
    %%%%%%%%%%%

    % Define the path to the 'sktparam' folder
    sktparam_folder = fullfile(outputFolder, 'sktparam');
    
    % Load the stored skew-t parameter objects
    load(fullfile(sktparam_folder, 'lc_skt.mat'), 'lc_skt');
    load(fullfile(sktparam_folder, 'sc_skt.mat'), 'sc_skt');
    load(fullfile(sktparam_folder, 'sh_skt.mat'), 'sh_skt');
    load(fullfile(sktparam_folder, 'df_skt.mat'), 'df_skt');
    
    % List of variable base names
    varnames = {'lc_skt', 'sc_skt', 'sh_skt', 'df_skt'};
    
    % Loop through and remove the covid+1 column if it exists
    for i = 1:length(varnames)
        varname = varnames{i};
        var = eval(varname);
    
        if size(var, 2) >= (covid + 1)
            var(:, covid + 1) = [];
        end
    
        % Save the updated variable back to the workspace
        assignin('base', varname, var);
    end
    
    % Generate a matrix of NaNs to store output
    zinf = NaN(size(actualvar));
    
    % check that pskt works
    % pskt(1, 0, 1, 0, 9999999) % cdf for normal distr at 1 0.8413
    
    % Loop through forecast dimensions
    for i = 1:size(lc_skt, 1)  % rows: horizons
        for j = 1:size(lc_skt, 2)  % columns: time
    
                % Evaluate pskt
                zinf(i, j) = pskt(actualvar(i, j), lc_skt(i, j), sc_skt(i, j), sh_skt(i, j), df_skt(i, j));
    
        end
    end

end


%% KS test

ksinf = zeros(size(zinf, 1), 1);
ksinfpv = zeros(size(zinf, 1), 1); 

for i = 1:size(zinf, 1) % fix the h-step ahead 

    valid_data = zinf(i, :)';   
    valid_data = valid_data(~isnan(valid_data));   % Remove missing values
     
    [ksinf(i), ksinfpv(i)] = kstestu(valid_data); % first compute the empirical CDF and then run the KS test 

end 

%% Rossi-Sekhposyan 2019

rvec = linspace(0, 1, 1000); % r is the support of the uniform distribution used to evaluate the PITs. 

% critical values are arranged 1, 5, 10% 
rs_boot_crit_value = repmat(struct('ks',[], 'cvm',[]), size(zinf,1), 1);

rs_test_stats = repmat(struct('ks_stat', [], 'cvm_stat', []), size(zinf,1), 1);
% if zero then implies test statistic below critical value at a given
% significance level => passes test. 
rs_logic = repmat(struct('ks_array', [], 'cvm_array', []), size(zinf,1), 1);

% not doing the p-values here since the distribution of the KS under RS
% 2019 is different from the regular KS, resulting in p-values unaligned
% with the rs_logic results. 


for i = 1:size(zinf, 1)

    z = zinf(i, :)';
    z = z(~isnan(z));

    [rs_boot_crit_value(i), rs_test_stats(i), rs_logic(i)] = rs_test(z, rvec);

end

% Currently the code exports the struct created above into a table. 
% Alternatively place each of the following things individually into the
% table

rs_ks_cv = vertcat(rs_boot_crit_value.ks); 
rs_cvm_cv = vertcat(rs_boot_crit_value.cvm);
rs_ks_test = vertcat(rs_test_stats.ks_stat);
rs_cvm_test = vertcat(rs_test_stats.cvm_stat);

% 0 = "accept" the null. Arranged by 1, 5, 10 cv. 
rs_ks_logic = vertcat(rs_logic.ks_array);
rs_cvm_logic = vertcat(rs_logic.cvm_array);


%% Berkowitz 2001

berK1inf= zeros(size(zinf, 1), 1); 
berK2inf= zeros(size(zinf, 1), 1); 
bert1inf= zeros(size(zinf, 1), 1);
bert2inf= zeros(size(zinf, 1), 1);

for i = 1:size(zinf, 1)

    [bert1inf(i), bert2inf(i), berK1inf(i), berK2inf(i) ] = berk(zinf(i, :)); % run Berkowitz with rho = 0 (ber1) and with \hat{rho} (ber2)

end 

%% Knueppel (2015)

% stet the parameters 
lags = -1;                % Let lag selection be automatic
prewhite = 0;             % No prewhitening

knueppel_stat= zeros(size(zinf, 1), 1); 
knueppel_pval= zeros(size(zinf, 1), 1); 

for i = 1:size(zinf, 1)
    
    [knueppel_stat(i), knueppel_pval(i)] = alpha0_1234_NW(zinf(i, :), lags, prewhite);

end 


%% Mitchell-Weale (2023)

% set the parameters (lags and prewhite are from Knueppel)
z_u = 0.95;         % Upper censoring threshold
z_l = 0.05;         % Lower censoring threshold
df = 5;             % Four moments and the frequency 

MW_stat= zeros(size(zinf, 1), 1); 
MW_pval= zeros(size(zinf, 1), 1); 

for i = 1:size(zinf, 1)
    
    stat_uncens = MW_alpha0_1234_NW(zinf(i, :), lags, prewhite, z_l, z_u); % same as knueppel but with ub and lb
    stat_freq = freqTestInCensoredRegion(zinf(i, :), lags, prewhite, z_l, z_u);
    MW_stat(i) = stat_uncens + stat_freq;
    MW_pval(i) = 1 - chi2cdf(MW_stat(i) , df);

end 

%% Save the data

% Define Excel filename once using the model name
filename = fullfile(evalDir, sprintf('%s_fcst_eval.xlsx', modellist{modelfcst + 1}));
horiz = (0:12)';

% Build table using the variables (assumed to be column vectors of length 13)
T = table(horiz, ...
          ksinf, ksinfpv, rs_ks_cv, rs_cvm_cv, ...
          rs_ks_test, rs_cvm_test, rs_ks_logic, rs_cvm_logic, ...
          bert1inf, bert2inf, ...
          berK1inf, berK2inf, ...
          knueppel_stat, knueppel_pval, ...
          MW_stat, MW_pval, ...
          'VariableNames', {'horiz', ...
                            'ksinf', 'ksinfpv', ...
                            'rs_ks_cv','rs_cvm_cv', ...
                            'rs_ks_test', 'rs_cvm_test', 'rs_ks_logic', 'rs_cvm_logic',...
                            'ber_1', 'ber_2', ...
                            'ber_1pv', 'ber_2pv', ...
                            'knueppel_stat', 'knueppel_pval', ...
                            'MW_stat', 'MW_pval'});

% Note that rs_ks_logic and rs_cvm_logic will produce FALSE and TRUE
% statements within the excel file. FALSE = fails to reject the null. 
% TRUE = accept the alternative.   

% Write table to Excel
writetable(T, filename);




