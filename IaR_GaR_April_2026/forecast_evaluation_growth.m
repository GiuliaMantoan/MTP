%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forecast evaluation code %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Authors: David Aikman, Rhys Bidder, Simon Lloyd, Guilia Mantoan, Simone
% Maso, Aditya Mori and Matthew Tong

close all; clear; clc;

%% Settings

set(0,'defaultAxesFontName', 'Times'); % font for chart axis
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1) % line style
rng('default'); % defaul random number generateor
rng(0);

cd 'C:\Users\k2368907\Dropbox\BoE-KCL Macro Forecasting\' % chage here the cd
addpath('Data\') % data input
addpath('Codes\intermediate_codes\')  % intemediate code (data loading..)
addpath('Codes\functions\')
addpath('Codes\functions\azzalini')  % skew t functions
outputFolder = fullfile(pwd, 'Outputs/');


% file names
fullFileName = fullfile('Data', 'GaRDataRaw_quarterly.xlsx'); % raw data

%% Import Data

%%%%%%%%%%%%%%%%%%%%%%
% Data loading panel %
%%%%%%%%%%%%%%%%%%%%%%


start_date = datetime('30-Sep-2007', 'InputFormat', 'dd-MMM-yyyy'); % first forecasted period
end_date   = datetime('31-Mar-2024', 'InputFormat', 'dd-MMM-yyyy'); % last date with actual data. 
end_fcst   = datetime('30-Mar-2027', 'InputFormat', 'dd-MMM-yyyy'); % forecast end date
covid_date = datetime('30-Jun-2020', 'InputFormat', 'dd-MMM-yyyy');
varnames = {'g4rgdp'}; % actual data variables you want to import. Here it is growth
ctrynames = {'UK'}; % actual data country you want to import. This will be the column from the sheet that is being imported.
momentlist = {'fcstmean', 'fcststdev', 'fcstskew'}; % moments for OLS and skewt

boefsctdata_growth; % this code focus only on growth

actualdata;

% actualvar %

%% Covid treatment

% drop the columns related to covid (all zeros) Q2-2020
covid = (year(covid_date) - year(start_date)) * 4 + (ceil(month(covid_date)/3) - ceil(month(start_date)/3));

if size(actualvar, 2) >= (covid +1)
    actualvar(:, (covid +1)) = [];
end

if size(mtestdata, 2) >= (covid*3 + 1)
    mtestdata(:, (covid*3 + 1):(covid*3 + 3)) = [];
end

%% Import the PITs

%%%%%%%%%%%%%%%%%%%
% Model selection %
%%%%%%%%%%%%%%%%%%%

modelfcst = 0; % 0 is boe, 1 is ols, 3 is skewt
modellist = {'boe', 'qreg'};


if modelfcst == 0

    %%%%%%%%%%%
    %   BOE   % (here we also compute mean, std dev and skew)
    %%%%%%%%%%%

    % Extract the corresponding columns from momgrowth (mode | mean | dispersion)
    modevar = mtestdata(:, 1:3:end);  % mode
    meanvar = mtestdata(:, 2:3:end); % mean
    square_root_dispersion_var = mtestdata(:, 3:3:end); % dispersion

    % generate a zero matrix to store the values
    zgrowth = NaN(size(actualvar));
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
                    zgrowth(i, j) = integral(@(y) ftp(y, m, s1, s2), -3, actualvar(i, j)); % get the cdf
                else
                    zgrowth(i, j) = NaN;
                end

                % get sd and skew
                stddevvar(i,j) = sqrt(sig);
                term1_skw = sqrt(2/pi) * (s2 - s1);
                term2_skw = (4/pi - 1) .* ((s2 - s1).^2 + s1 .* s2);
                skewvar(i,j) =  term1_skw .* term2_skw;

            else

                if ~isnan(actualvar(i, j))
                    zgrowth(i, j) = normcdf((actualvar(i, j) - meanvar(i,j)) / square_root_dispersion_var(i,j)); % standardized value for growth if mean = mode
                else
                    zgrowth(i, j) = NaN;
                end

                % get sd and skew
                stddevvar(i,j) = square_root_dispersion_var(i,j);
                skewvar(i,j) =  0;

            end
            % end
        end
    end

    % make sure no NaN col
    zgrowth = zgrowth(:, ~all(isnan(zgrowth), 1));

    % save the outcome for Boe Forecast
    fcstmean   = [meanvar(:,1:covid),   zeros(size(meanvar,1),1),   meanvar(:,covid+1:end)];
    fcststdev = [stddevvar(:,1:covid),  zeros(size(stddevvar,1),1), stddevvar(:,covid+1:end)];
    fcstskew   = [skewvar(:,1:covid),   zeros(size(skewvar,1),1),   skewvar(:,covid+1:end)];

    %% ADDED !! 
    % The below section has been added because previously there was an
    % issue in line 160 whereby currentVar featured 52 columns despite
    % there being only 51. This arises because of lines 136-139 which adds
    % covid+1. The snippet below strips the "+1". 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fcstmean  = fcstmean(:, 1:covid);
    % fcststdev = fcststdev(:,1:covid);
    % fcstskew  = fcstskew(:, 1:covid);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % settings
    quarterDates = start_date:calmonths(3):end_date;
    momentlist = {'fcstmean', 'fcststdev', 'fcstskew'};
    horizons = 13;
    filename = fullfile(outputFolder, 'boe_growth.xlsx');

    % Loop over each variable in the list
    for i = 1:length(momentlist)

        % Get the current forecast data
        currentVar = eval(momentlist{i});

        % Create the table with quarterly dates as the first column
        T = [table(quarterDates', 'VariableNames', {'Dates'}), array2table(currentVar')];

        % Rename the forecast columns from the second column onward as h_0, h_1, ..., h_12
        forecastNames = strcat("h_", string(0:horizons-1));
        T.Properties.VariableNames(2:end) = cellstr(forecastNames);
        % Export the table to a specific sheet in the Excel file
        writetable(T, filename, 'Sheet', momentlist{i});

    end

    clearvars stddevvar  term1_skw  term2_skw  skewvar  fcstmean fcststdev fcstskew quarterDates  momentlist horizons filename currentVar T forecastNames

elseif modelfcst == 1
    
    

end

%% KS test

ksgrowth = zeros(size(zgrowth, 1), 1);
ksgrowthpv = zeros(size(zgrowth, 1), 1);

for i = 1:size(zgrowth, 1) % fix the h-step ahead

    valid_data = zgrowth(i, :)';
    valid_data = valid_data(~isnan(valid_data));   % Remove missing values

    [ksgrowth(i), ksgrowthpv(i)] = kstestu(valid_data); % first compute the empirical CDF and then run the KS test

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
KS = size(zinf, 1);
CVM = size(zinf, 1);

for i = 1:size(zinf, 1)

    z = zinf(i, :)';
    z = z(~isnan(z));

    [rs_boot_crit_value(i), rs_test_stats(i), rs_logic(i), KS(i), CVM(i)] = rs_test(z, rvec);

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

berK1growth= zeros(size(zgrowth, 1), 1);
berK2growth= zeros(size(zgrowth, 1), 1);
bert1growth= zeros(size(zgrowth, 1), 1);
bert2growth= zeros(size(zgrowth, 1), 1);

for i = 1:size(zgrowth, 1)

    [bert1growth(i), bert2growth(i), berK1growth(i), berK2growth(i) ] = berk(zgrowth(i, :)); % run Berkowitz with rho = 0 (ber1) and with \hat{rho} (ber2)

end

%% Knueppel (2015)

% stet the parameters
lags = -1;                % Let lag selection be automatic
prewhite = 0;             % No prewhitening

knueppel_stat= zeros(size(zgrowth, 1), 1);
knueppel_pval= zeros(size(zgrowth, 1), 1);

for i = 1:size(zgrowth, 1)

    [knueppel_stat(i), knueppel_pval(i)] = alpha0_1234_NW(zgrowth(i, :), lags, prewhite);

end


%% Mitchell-Weale (2023)

% set the parameters (lags and prewhite are from Knueppel)
z_u = 0.95;         % Upper censoring threshold
z_l = 0.05;         % Lower censoring threshold
df = 5;             % Four moments and the frequency

MW_stat= zeros(size(zgrowth, 1), 1);
MW_pval= zeros(size(zgrowth, 1), 1);

for i = 1:size(zgrowth, 1)

    stat_uncens = MW_alpha0_1234_NW(zgrowth(i, :), lags, prewhite, z_l, z_u); % same as knueppel but with ub and lb
    stat_freq = freqTestInCensoredRegion(zgrowth(i, :), lags, prewhite, z_l, z_u);
    MW_stat(i) = stat_uncens + stat_freq;
    MW_pval(i) = 1 - chi2cdf(MW_stat(i) , df);

end

%% Galvao, Mantoan and Mitchell

z = zgrowth';

MC = 1000; %Number of MC replication for power exercise
bootMC = 1000;%Number of MC replication for bootstrap CV
randn('seed',bootMC); %#ok<*RAND>
rand('seed',bootMC);  %#ok<*RAND>
randn(MC,1);
rand(MC,1);

tic
rvec = 0:0.001:1;

% Number of horizons
%H=[4, 8, 12, 24, 50, 100];
%H = [5, 9, 13]; % because nowcast is horizon 0.
talesstackks = nan(1,3);
talesstackcvm = nan(1,3);
tableks_bonf= nan(1,3);
tablecvm_bonf= nan(1,3);

P = size(z, 1);
%R = Rinvec;
%T = R + P;

% loop over the different horizons
%for h=1:size(H,2)
% Hsel=H(h);
% disp(Hsel)

el = floor (P^(1/4));

QVrejvecs= zeros(MC,3);
CVMrejvecs= zeros(MC,3);
QVrejvecs_bonf= zeros(MC,3);
CVMrejvecs_bonf= zeros(MC,3);
Hz=size(z,2);

parfor j = 1:MC
    stream1 = RandStream('mrg32k3a','seed',4829575);
    stream1.Substream = j;
    [QVrejvecsj,CVMrejvecsj,QVrejvecs_bonfj,CVMrejvecs_bonfj]=size_statistic_h2(z, KS, CVM, Hz ,stream1,rvec,el, bootMC) ;
    QVrejvecs(j,:)=QVrejvecsj;
    CVMrejvecs(j,:)=CVMrejvecsj;
    %QVrejvecs_h(j,:)=QVrejvecs_hj;
    %CVMrejvecs_h(j,:)=CVMrejvecs_hj;
    QVrejvecs_bonf(j,:)=QVrejvecs_bonfj;
    CVMrejvecs_bonf(j,:)=CVMrejvecs_bonfj;
    %QVrejvecs_rs(j,:)=QVrejvecsj_rs;
    %CVMrejvecs_rs(j,:)=CVMrejvecsj_rs;
end

%tableksh=struct(append('h','_',num2str(h)),mean(QVrejvecs_h,1));
%tablecvmh=struct(append('h','_',num2str(h)),mean(CVMrejvecs_h,1));
tableks_bonf=mean(QVrejvecs_bonf,1);
tablecvm_bonf=mean(CVMrejvecs_bonf,1);
talesstackks=mean(QVrejvecs,1);
talesstackcvm=mean(CVMrejvecs,1);
%talesks_rs(h,:)=mean(QVrejvecs_rs,1);
%talescvm_rs(h,:)=mean(CVMrejvecs_rs,1);

timeElapsed = toc;
%end

tab_maxsup = [
    talesstackks(:,1),  tableks_bonf(:,1), ...
    talesstackks(:,2),  tableks_bonf(:,2), ...
    talesstackks(:,3),  tableks_bonf(:,3);
    talesstackcvm(:,1), tablecvm_bonf(:,1), ...
    talesstackcvm(:,2), tablecvm_bonf(:,2), ...
    talesstackcvm(:,3), tablecvm_bonf(:,3)
    ];

% rows 1-3: KS test statistic at horizons 4, 8 and 12
% rows 4-6: CvM test statistic at horizons 4, 8 and 12
% columns 1, 4, 7: standard sup tests with three weights - none, w, and
% alternative wi
% columns 2, 5, 8: robust sup tests with the three weighting schemes
% columns 3, 6, 9: Bonferroni-adjusted tests with the three weightings.

% each of the entries tells you frequency of rejects so tab_maxsup(1,1) =
% KS test stat is rejected 48% of the time.

tab_maxsup2= [talesstackks(:,1), tableks_bonf(:,1) talesstackcvm(:,1)  tablecvm_bonf(:,1);
    talesstackks(:,2)  tableks_bonf(:,2) talesstackcvm(:,2) tablecvm_bonf(:,2);
    talesstackks(:,3)  tableks_bonf(:,3) talesstackcvm(:,3)  tablecvm_bonf(:,3)];


%% Save the data

% Define Excel filename once using the model name
filename = fullfile(outputFolder, sprintf('%s_fcst_eval_growth.xlsx', modellist{modelfcst + 1}));
horiz = (0:12)';

% Build table using the variables (assumed to be column vectors of length 13)
T = table(horiz, ...
    ksgrowth, ksgrowthpv, rs_ks_cv, rs_cvm_cv, ...
    rs_ks_test, rs_cvm_test, rs_ks_logic, rs_cvm_logic, ...
    bert1growth, bert2growth, ...
    berK1growth, berK2growth, ...
    knueppel_stat, knueppel_pval, ...
    MW_stat, MW_pval, ...
    'VariableNames', {'horiz', ...
    'ksgrowth', 'ksgrowthpv', ...
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

rowNames = {'KS','CVM'};
varNames = { 'Standard', 'Bonf', 'Standard_w', 'Bonferroni_w', 'Standard_wi', 'Bonferroni_wi'};

% saving the GMM test results
gmmFilename = fullfile(outputFolder, sprintf('GMM_growth_%s.xlsx', modellist{modelfcst + 1}));

GMM_T    = array2table(tab_maxsup, 'RowNames', rowNames, 'VariableNames', varNames);

writetable(GMM_T, gmmFilename, 'Sheet', 'MaxSup', 'WriteRowNames', true);

% write tab_maxsup2_tbl to sheet “MaxSup2”
writematrix(tab_maxsup2, gmmFilename, 'Sheet', 'MaxSup2');



