% clean
close all; clear; clc;

opts = spreadsheetImportOptions("NumVariables", 6);
opts.Sheet = "UK_monthly";
opts.DataRange = "A2:F666";
opts.VariableNames = ["ccode", "date", "g4Infl", "Oil", "ue", "vu_tight"];
opts.VariableTypes = ["categorical", "string", "double", "double", "double", "double"];
opts = setvaropts(opts, "date", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ccode", "date"], "EmptyFieldRule", "auto");
kf_data = readtable("C:\Users\simon\Dropbox\BoE-KCL Macro Forecasting\Data\others\NAIRU_est.xlsx", opts, "UseExcel", false);
clear opts

addpath('C:\Users\simon\Dropbox\BoE-KCL Macro Forecasting\Codes\functions\Kalman_filter')

%% Converting data format.
% Data is in strings
months = kf_data.date;
% Extract year and quarter number
year = str2double(extractBefore(months, 'm'));
months = str2double(extractAfter(months, 'm'));
% Construct datetime from year and month (use 1st day of month)
kf_data.date_dt = datetime(year, months, 1);

%% HP Filter

%{

% Identify the last index in Q4 1996
start_date    = datetime(1971,2,1);  % first obs available
est_date      = datetime(1988,1,1);  % first point at which the filtered series is produced
end_date      = datetime(2025,3,1);  % last obs available

idxstart       = find(kf_data.date_dt <= start_date, 1, 'last'); % idx of first obs
idxest       = find(kf_data.date_dt <= est_date, 1, 'last'); % idx of first point at which filtered series is produced
N       = find(kf_data.date_dt <= end_date, 1, 'last'); % last obs of the var I want to filter

series_list = {'ue', 'vu_tight'};

for k = 1:numel(series_list)
    
    % get series name
    series_name = series_list{k};

    % extract the series
    y = kf_data.(series_name);  
    
    % for gdp take the log 
    %if k == 2
    %   y = log(y);
    %end

    % size of (iterated) HP filtered series
    nCols = N - idxest + 1;

    % Pre-allocate the N×T matrix
    hp_one_sided = NaN(N, nCols);

    lambda = 129600;  % smoothing parameter (originally we use 14400 - then following Uhlig-Ravn we set it to 
    
    for j = 1 : nCols

        t = idxest + j - 1;       % index of this terminal quarter
        ysub = y(idxstart:t);             % subsample up to t
    
        % two-sided hpfilter on subsample; at the end-point this
        % effectively uses only past data
        [trend_sub, ~] = hpfilter(ysub,"Smoothing", lambda);
    
        % store the endpoint trend in row t, column j
        hp_one_sided(t, j) = trend_sub(end);

    end

    % acess the filtered series
    trend_series(:, k) = diag( hp_one_sided(idxest:end, :) ); % for Simone

    % gap 
    gap_series(:,k) = y(idxest:N) - trend_series(:, k); % copy paste inside GaRData

end

% checks sm 
% series at which you you want to apply the HP filter
%series_list = {'ue'};
%Y = kf_data.ue(idxstart:idxest);
%[trend_sub, ~] = hpfilter(Y, 14400);

%}


%% 3.  Extracting the series a
infl   = kf_data.g4Infl;                         % CPI inflation
oil    = kf_data.Oil;                            % oil price
u      = kf_data.ue;                             % unemployment rate (%)

valid = ~isnan(u);

infl = infl(valid);
oil  = oil(valid);
u    = u(valid);

% first‑differences
Dpi       = [NaN; diff(infl)];
oilLog    = log(oil);
oilInfl   = [NaN; diff(oilLog)];                 % NaN at t=1, then %∆ oil

%Dz     = [NaN; diff(oil)];
Dz      = [NaN; diff(oilInfl)];

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

valid1   = all(~isnan([Dpi, drivers1]),2);  
Y1       = Dpi(valid1);
X1       = drivers1(valid1,:);
dd = kf_data.date_dt(valid);
Dates1   = dd(valid1);
T1       = numel(Y1);

%% Set up for Specification 2
drivers2 = [Dpi_l1, Dpi_l2, Dz_l0, Dz_l1, Dz_l2, U_l0]; 

valid2   = all(~isnan([Dpi, drivers2]),2);
Y2       = Dpi(valid2);
X2       = drivers2(valid2,:);
dd       = kf_data.date_dt(valid);
Dates2   = dd(valid2);
T2       = numel(Y2);

%%  Maximum Likelihood Estimation
varianceRatio = 0.16;          % var(η) / var(ε) - arbitary restriction of signal to noise ratio.

% Specification 1
theta0_1 = [0.5 0.2 0.30 0.15 0.10 0.05  log(0.3^2)  mean(X1(:,5))]; % taking logs to ensure that the variance is non-negative
%           a1  a2   b1   g1   g2    vareps      u
logLik1  = @(th) kalmanPhillips_lag(th,Y1,X1,varianceRatio);

% could do SQP to have constrained optimisation but corner solutions
% sometimes arise
optsUNC  = optimoptions('fmincon','algorithm','active-set','display','iter', ...
    'maxiterations',400,'TolFun',1e-8,'TolX',1e-8);

testval = logLik1(theta0_1);
disp(['Initial log-likelihood = ', num2str(testval)]);

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
theta0_2 = [0.5 0.2 0.30 0.15 0.10 0.05  log(0.3^2)  mean(X2(:,6))];
logLik2  = @(th) kalmanPhillips_contemp(th,Y2,X2,varianceRatio);

testval = logLik2(theta0_2);
disp(['Initial log-likelihood = ', num2str(testval)]);

[theta2_hat,nLL2,flag2] = fminunc(logLik2,theta0_2,optsUNC);

% unpack
[a1_2,a2_2,beta2,g20,g21,g22,lvar_eps2,x0_2] = deal(theta2_hat(1),theta2_hat(2), ...
    theta2_hat(3),theta2_hat(4), ...
    theta2_hat(5),theta2_hat(6), ...
    theta2_hat(7),theta2_hat(8));

var_eps2 = exp(lvar_eps2);   
var_eta2 = varianceRatio*var_eps2;

[~,~,xf2,~,xs2,~] = kalmanPhillips_contemp(theta2_hat,Y2,X2,varianceRatio);

nairu_spec1 = xf1;
nairu_spec2 = xf2;
xf1_lag_removed = xf1(2:end);
nairu_spec_avg = (xf2 + xf2)/2; % for Simone

nairu_rts1 = xs1;
nairu_rts2 = xs2;
xs1_lag_removed = xs1(2:end);
nairu_rts_avg = (xs2 + xs2)/2;

% construct unemployment gap (note that valid1 = valid2)
ugap_kf = u(valid1) - nairu_spec_avg;

%% Result Display
fprintf('\n===========  Specification 1  (Lagged Gap)  ===========\n');
fprintf(' a1=%6.3f  a2=%6.3f  β=%6.3f  γ1=%6.3f  γ2=%6.3f\n',a1_1,a2_1,beta1,g11,g12);
fprintf(' σ_ε=%6.3f  σ_η=%6.3f   logL=%.2f  (flag %d)\n',sqrt(var_eps1),sqrt(var_eta1),-nLL1,flag1);

fprintf('\n===========  Specification 2  (Contemporaneous Gap)  ===========\n');
fprintf(' a1=%6.3f  a2=%6.3f  β=%6.3f  γ0=%6.3f  γ1=%6.3f  γ2=%6.3f\n', ...
    a1_2,a2_2,beta2,g20,g21,g22);
fprintf(' σ_ε=%6.3f  σ_η=%6.3f   logL=%.2f  (flag %d)\n',sqrt(var_eps2),sqrt(var_eta2),-nLL2,flag2);


%% Plots
figure('Name','NAIRU');
plot(Dates1,u(valid1),'k','LineWidth',1); hold on;
plot(Dates1,nairu_spec_avg,'r','LineWidth',1.5);
plot(Dates2,nairu_rts_avg,'b','LineWidth',1.5);
legend({'Observed unemployment','NAIRU Specification Average','NAIRU Smoothed Average'},'Location','best');
ylabel('Unemployment rate (%)'); grid on;

% notes:
%Warning: You have passed FMINCON options to FMINUNC. FMINUNC will use the common options and ignore the FMINCON options that do not
%apply.
% That warning means you built an options object for the constrained nonlinear solver (fmincon)
% and then handed it to the unconstrained solver (fminunc). fminunc only understands its own subset of options