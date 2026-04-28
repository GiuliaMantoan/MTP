                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %%% Extra series construction %%%
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Authors: David Aikman, Rhys Bidder, Simone Maso, Aditya Mori 
% First version: May 2025
% This version: May 2025

% NOTE: this code perform the following operations:
% 1) extract the trend from HP filtering the output series 
% 2) extract the trend from HP filtering the unemployment series 
% 3) get the filtered series for unemployment using kalman filter
% 4) compute the exponential weighted moving average inflation 

%% Set up the Import Options and import the data

% clean
close all; clear; clc;

% note: with resepct to the original dataset provided by BOE
% GaRDataRaw.xlxs (sheet: MAIN data) that ends in 2018Q4, we extend it to
% 2022Q4 for UK

opts = spreadsheetImportOptions("NumVariables", 24);

% Specify sheet and range
opts.Sheet = "MAIN data";
opts.DataRange = "A2:X3213";

% Specify column names and types
opts.VariableNames = ["ccode", "date", "rgdp_oecd", "tce_linear", "dep", "idfci_nocredhp", "gr3y_realhp", "diff3y_creditgdp", "ca_smooth", "inflation", "diff1y_cbrate", "volatility", "ue", "gr3y_realtcredit", "diff1y_creditgdp", "gr3y_tcredit", "hh_3y", "PNFC_3y", "creditgdp", "gr1y_realhp", "Bank_3y", "Nonbank_3y", "rgdp", "g4Infl"];
opts.VariableTypes = ["categorical", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "date", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ccode", "date"], "EmptyFieldRule", "auto");

% Import the data
data = readtable("C:\Users\k2370179\Dropbox\BoE-KCL Macro Forecasting\Data\GaRDataRaw.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts

%% Converting date format 

% Assuming your variable is a string array
qtrs = data.date;  % e.g., ["1972q1", "1972q2", ...]

% Extract year and quarter number
year = str2double(extractBefore(qtrs, 'q'));
quarter = str2double(extractAfter(qtrs, 'q'));

% Convert quarter to starting month (Q1 = Jan, Q2 = Apr, Q3 = Jul, Q4 = Oct)
month = (quarter - 1) * 3 + 1;

% Construct datetime from year and month (use 1st day of month)
data.date_dt = datetime(year, month, 1);

% Keep only UK data
uk_data = data(data.ccode=='UK',:);

% Identify the last index in Q4 1996
N = height(uk_data);
cutoffDate      = datetime(1990,1,1);         % CHANGE HERE DEPENDING ON THE CUT-OFF DATE
idxCutoff       = find(uk_data.date_dt <= cutoffDate, 1, 'last');

% series at which you you want to apply the HP filter
series_list = {'ue', 'rgdp'};

%% HP Filter

for k = 1:numel(series_list)
    
    % get series name
    series_name = series_list{k};

    % extract the series
    y = uk_data.(series_name);  
    
    % for gdp take the log 
    if k == 2
       y = log(y);
    end

    % start cut-off
    startIdx = idxCutoff; 
    nCols = N - startIdx + 1;

    % Pre-allocate the N×T matrix
    hp_one_sided = NaN(N, nCols);

    lambda = 1600;  % smoothing parameter
    
    for j = 1 : nCols
        t = startIdx + j - 1;       % index of this terminal quarter
        ysub = y(1:t);             % subsample up to t
    
        % two-sided hpfilter on subsample; at the end-point this
        % effectively uses only past data
        [trend_sub, ~] = hpfilter(ysub, lambda);
    
        % store the endpoint trend in row t, column j
        hp_one_sided(t, j) = trend_sub(end);
    end

    % This starts from Q4 1996 and no time period before. 
    trend_series(:, k) = diag( hp_one_sided(startIdx:end, :) ); % for Simone

    % gap 
    gap_series(:,k) = y(startIdx:end) - trend_series(:, k); % copy paste inside GaRData

end

%% Exponentially Weighted Average 

% get inflation
infl = uk_data.g4Infl;  

% pre-allocation and initializiation
ewa_inf = NaN(size(infl));
ewa_inf(1) = infl(1);

% smoothign
alpha = 0.25;

% Recursive computation
for t = 2:N
    ewa_inf(t) = alpha * infl(t) + (1 - alpha) * ewa_inf(t-1);
end


%% Set up the Import Options and import the data

% clean
close all; clear; clc;

opts = spreadsheetImportOptions("NumVariables", 25);
opts.Sheet = "UK";
opts.DataRange = "A2:Y205";
opts.VariableNames = ["ccode", "date", "rgdp_oecd", "tce_linear", "dep", "idfci_nocredhp", "gr3y_realhp", "diff3y_creditgdp", "ca_smooth", "inflation", "diff1y_cbrate", "volatility", "ue", "gr3y_realtcredit", "diff1y_creditgdp", "gr3y_tcredit", "hh_3y", "PNFC_3y", "creditgdp", "gr1y_realhp", "Bank_3y", "Nonbank_3y", "oil", "exports", "g4Infl"];
opts.VariableTypes = ["categorical", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, "date", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ccode", "date"], "EmptyFieldRule", "auto");
kf_data = readtable("C:\Users\k2370179\Dropbox\BoE-KCL Macro Forecasting\Data\others\NAIRU_est.xlsx", opts, "UseExcel", false);
clear opts

%% Converting data format.
% Data is in strings
qtrs = kf_data.date;
% Extract year and quarter number
year = str2double(extractBefore(qtrs, 'q'));
quarter = str2double(extractAfter(qtrs, 'q'));
% Convert quarter to starting month (Q1 = Jan, Q2 = Apr, Q3 = Jul, Q4 = Oct)
month = (quarter - 1) * 3 + 1;
% Construct datetime from year and month (use 1st day of month)
kf_data.date_dt = datetime(year, month, 1);

%% 3.  Extracting the series a
infl   = kf_data.g4Infl;                      % CPI inflation
oil    = kf_data.oil;                            % oil price
u      = kf_data.ue;                             % unemployment rate (%)

% first‑differences
Dpi    = [NaN; diff(infl)];
oilLog    = log(oil);
oilInfl   = [NaN; diff(oilLog)];                 % NaN at t=1, then %∆ oil

%Dz     = [NaN; diff(oil)];
Dz     = [NaN; diff(oilInfl)];

% lags
% inflation
Dpi_l1 = lagmatrix(Dpi,1);
Dpi_l2 = lagmatrix(Dpi,2);
% oil prices
Dz_l0  = Dz;                 % contemporaneous (for Spec 2)
Dz_l1  = lagmatrix(Dz,1);
Dz_l2  = lagmatrix(Dz,2);
% unemployment
U_l0   = u;                  % contemporaneous (Spec 2)
U_l1   = lagmatrix(u,1);     % lagged         (Spec 1)

%% Set up for Specification 1
drivers1 = [Dpi_l1, Dpi_l2, Dz_l1, Dz_l2, U_l1];     % 5 columns

valid1   = all(~isnan([Dpi, drivers1]),2);           % restricts sample size to be the same across all variables.
Y1       = Dpi(valid1);
X1       = drivers1(valid1,:);                       % matrix and every element of it is taken
Dates1   = kf_data.date_dt(valid1);
T1       = numel(Y1);

%% Set up for Specification 2
drivers2 = [Dpi_l1, Dpi_l2, Dz_l0, Dz_l1, Dz_l2, U_l0]; % 6 columns (γ0,γ1,γ2)

valid2   = all(~isnan([Dpi, drivers2]),2);
Y2       = Dpi(valid2);
X2       = drivers2(valid2,:);
Dates2   = kf_data.date_dt(valid2);
T2       = numel(Y2);

%%  Maximum Likelihood Estimation
varianceRatio = 0.16;          % var(η) / var(ε) - arbitary restriction of signal to noise ratio.

% Specification 1
theta0_1 = [0.5 0.2 0.30 0.10 0.05  log(0.3^2)  mean(u)]; % taking logs to ensure that the variance is non-negative
%           a1  a2   b1   g1   g2    vareps      u*
logLik1  = @(th) kalmanPhillips_lag(th,Y1,X1,varianceRatio);

% could do SQP to have constrained optimisation but corner solutions
% sometimes arise
optsUNC  = optimoptions('fminunc','algorithm','quasi-newton','display','iter', ...
    'maxiterations',400,'TolFun',1e-8,'TolX',1e-8);

[theta1_hat,nLL1,flag1] = fminunc(logLik1,theta0_1,optsUNC);

% unpack
[a1_1,a2_1,beta1,g11,g12,lvar_eps1,x0_1] = deal(theta1_hat(1),theta1_hat(2), ...
    theta1_hat(3),theta1_hat(4), ...
    theta1_hat(5),theta1_hat(6), ...
    theta1_hat(7));

var_eps1 = exp(lvar_eps1);
var_eta1 = varianceRatio*var_eps1;

% filter/smooth with final θ̂
[~,~,xf1,~,xs1,~] = kalmanPhillips_lag(theta1_hat,Y1,X1,varianceRatio);

% Specification 2
theta0_2 = [0.5 0.2 0.30 0.15 0.10 0.05  log(0.3^2)  mean(u)];
logLik2  = @(th) kalmanPhillips_contemp(th,Y2,X2,varianceRatio);
[theta2_hat,nLL2,flag2] = fminunc(logLik2,theta0_2,optsUNC);

% unpack
[a1_2,a2_2,beta2,g20,g21,g22,lvar_eps2,x0_2] = deal(theta2_hat(1),theta2_hat(2), ...
    theta2_hat(3),theta2_hat(4), ...
    theta2_hat(5),theta2_hat(6), ...
    theta2_hat(7),theta2_hat(8));
var_eps2 = exp(lvar_eps2);   var_eta2 = varianceRatio*var_eps2;

[~,~,xf2,~,xs2,~] = kalmanPhillips_contemp(theta2_hat,Y2,X2,varianceRatio);

nairu_spec1 = xf1;
nairu_spec2 = xf2;
nairu_spec_avg = (xf1 + xf2)/2; % for Simone

nairu_rts1 = xs1;
nairu_rts2 = xs2;
nairu_rts_avg = (xs1 + xs2)/2;

% construct unemployment gap (note that valid1 = valid2)
ugap_kf = u(valid1) - nairu_spec_avg;

%% Plots
%figure('Name','NAIRU');
%plot(Dates1,u(valid1),'k','LineWidth',1); hold on;
%plot(Dates1,nairu_spec_avg,'r','LineWidth',1.5);
%plot(Dates2,nairu_rts_avg,'b','LineWidth',1.5);
%legend({'Observed unemployment','NAIRU Specification Average','NAIRU Smoothed Average'},'Location','best');
%ylabel('Unemployment rate (%)'); grid on;

