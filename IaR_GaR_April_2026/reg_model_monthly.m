                                    %%%%%%%%%%%%%%%%%%%%%%%%%
                                    %%% Regression models %%%
                                    %%%%%%%%%%%%%%%%%%%%%%%%%

% Authors: David Aikman, Rhys Bidder, Simone Maso, Aditya Mori 
% First version: May 2025
% This version: July 2025

close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS AND PATHS %
%%%%%%%%%%%%%%%%%%%%%%%

set(0,'defaultAxesFontName', 'Times'); % font for chart axis 
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1) % line style 
rng('default'); % defaul random number generateor 

cd 'C:\Users\simon\Dropbox\BoE-KCL Macro Forecasting\' % change here the cd 
%cd 'C:\Users\336770\Dropbox\BoE-KCL Macro Forecasting\' % change here the cd 
addpath('Data\') % data input
addpath('Codes\intermediate_codes\')  % intermediate code (data loading..)
addpath('Codes\functions\')   % functions
addpath('Codes\functions\azzalini')  % skew t functions 
addpath('Codes\functions\CRPS')  %  Continuos rank prob score functions 
addpath('Codes\functions\Simon_qreg')  %  Continuos rank prob score functions 
addpath('Codes\functions\two_piece_normal')  %  TPN functions
outputFolder = fullfile(pwd, 'Outputs/'); % charts and tables
dropboxFolder = 'C:\Users\simon\Dropbox\Apps\Overleaf\BoE-KCL Macro Forecasting\TablesFigures'; % chart directory (SM dropbox) 
%dropboxFolder = 'C:\Users\336770\Dropbox\Apps\Overleaf\BoE-KCL Macro Forecasting\TablesFigures'; % chart directory (GM dropbox) 

% File names
fullFileName = fullfile('Data', 'IaRDataRaw_monthly.xlsx'); % raw data 

%%%%%%%%%%%%%%%%%%
%% CONTROL PANEL %
%%%%%%%%%%%%%%%%%%

startT = datenum(1990,02,01);  %  --> 1990m1
endT   = datenum(2025,01,01);  %  --> 2024m12

ctrynames = {'UK'}; % country 
onlyuk      = 1; % Set to 1 for UK only, 0 for panel (multiple countries to be specified in ctrynames) [For now only UK]

dep_var_name = {'g4cpi'}; % dep var inflation (g4Infl) - for CPI is (g4cpi) - for CPI ex energy is (g4cpinoener)
infl_persistence_cat = {'avgcpi'}; % inflation persistence (avg12Infl) (average inflation over the past year) - for CPI is (avgcpi) - for CPI ex energy is (avgcpinoener)
infl_expect_cat = {'infl1_m'}; % Inflation expectations (next calendar year average inflation expectation)
economic_slack_cat = {'delta1vu_tightness'}; % Economic slack measure 1) delta 1y labour market tightness: delta1vu_tightness  2) unemployment gap hp filter (ugap_hp_filter_lambda_129600) 3) unemployment gap (ugap_kalman_filter) 4) v/u tightness (vu_tightness)
supply_shock_cat = {'yoy_growth_import_defl_FAME'}; % supply shock measure 1) Oil price inflation (g4oil) 2) G7 core inflation (g4cpicore_global) 3) Import deflator (yoy_growth_import_defl_FAME)
fci_cat = {'bond_spread'}; % financial conditions  
% economic_slack_level_cat = {'vu_tightness'}; 

model_selection = 4; % 1 is OLS + normal, 2 is quantile + skewt, 3 quantile + semiparam density (Mitchell, Poon, Zhu), 4 is quantile + two piece normal  [For now only model 4 available]
modellist = {'ols', 'skewt','semi-param', 'two-piece-normal'}; % name of the models [for now we do not use model 1: OLS]
horizons = 37; % forecast horizons
momentlist = {'fcstmean', 'fcststdev', 'fcstskew'}; % moments 
quantilelevels = 0.05:0.05:0.95; % quantiles for qreg 
StartEst =  datenum(2004,01,01); % first period at which fcst are estimated datenum(2004,01,01)
dummyvarname = {'covid'};
covid = 0; % 1 covid is accounted (dummy = 1 for 2020q2 and 2020q3) 0 covid not accounted 
covid_dates = datenum(2020,04,01); % first date at which covid appears (if you want to change the full structure of covid dummies you need to go in GARdataRaw.xlsx (sheet 'covid')

% bootstrap option Simon
bst = 1;
bstOptions.blocksize = 24; 
bstOptions.nboot    = 10; % Usually 5000
bstOptions.ci       = 68;  % just for the code to run - it never get used

horizons_crps = [1 13 25]; % horizons of which you want to know the best spec in terms of crps (current | 1y | 2y)
quantilesirf = [0.1 0.5 0.9]; % quantiles to plot in the IRF chart
use_best_spec = 1; % plot IRF for the best specification (in stats terms according to the CRPS) within a given model for the last time period available + save the moments (mean, sd and skew) for the best spec
specplot = 0; % if use_best_spec = 0 and you want to see the IRF and save the moments of a given specification just select here the spec of interest within a given model
plot_bands_irf = 1; % 1 68% bands; 0 no bands

%%%%%%%%%%%%%%%%%%%%%%
%% ACTUAL OPERATIONS %
%%%%%%%%%%%%%%%%%%%%%%

%% Generate the different model combinations

%[i0,i1,i2,i3,i4,i5,i6] = ndgrid(1:size(dep_var_name,2), 1:size(infl_persistence_cat,2), 1:size(infl_expect_cat,2), 1:size(economic_slack_cat,2), 1:size(supply_shock_cat,2), 1:size(fci_cat,2), 1:size(economic_slack_level_cat,2));
[i0,i1,i2,i3,i4,i5] = ndgrid(1:size(dep_var_name,2), 1:size(infl_persistence_cat,2), 1:size(infl_expect_cat,2), 1:size(economic_slack_cat,2), 1:size(supply_shock_cat,2), 1:size(fci_cat,2));

% Put the variables in array
vars = {dep_var_name(i0(:)), infl_persistence_cat(i1(:)), infl_expect_cat(i2(:)), economic_slack_cat(i3(:)), supply_shock_cat(i4(:)), fci_cat(i5(:))};

% Loop over each and transpose if it has more than one column
for k = 1:numel(vars)
    v = vars{k};
    if size(v, 2) > 1
        vars{k} = v.';
    end
end

% get the combos
combo_specifications = [vars{1}, vars{2}, vars{3}, vars{4}, vars{5}, vars{6}];

%% pre-allocate the results 

dv1 = datevec(StartEst); dv2 = datevec(endT); % how many months between start est and end date 
nMonths = (dv2(1) - dv1(1))*12 + (dv2(2) - dv1(2)) + 1;

coeffqr_SL = zeros(size(vars,2) +1,size(quantilelevels,2),horizons,nMonths,size(combo_specifications,1));
bootstrapqrg_SL = zeros(size(vars,2) +1,size(quantilelevels,2),horizons, bstOptions.nboot,size(combo_specifications,1));
        
% empty vector to store the predicted quantile 
predicted_quantiles_SL = zeros(nMonths, size(quantilelevels,2), horizons,size(combo_specifications,1));

% semi parametric distribution
semi_param_distr = NaN(nMonths, 20000, horizons, size(combo_specifications,1));
empirical_cdf_semi_param = NaN(20000+1, 2, nMonths, horizons, size(combo_specifications,1)); % dim 2 is 1 for the support and the other for the cdf

% skew t distribution
lc_skt = NaN(nMonths,horizons,size(combo_specifications,1));
sc_skt = NaN(nMonths,horizons,size(combo_specifications,1)); 
sh_skt = NaN(nMonths,horizons,size(combo_specifications,1)); 
df_skt = NaN(nMonths,horizons,size(combo_specifications,1));
% param_skt = NaN(nMonths,4,horizons,size(combo_specifications,1));

% two piece normal 
param_tpn = NaN(nMonths, 3, horizons,size(combo_specifications,1)); % mode, s1 and s2

% moments - for all distributions
fcstmean = NaN(nMonths,horizons,size(combo_specifications,1));
fcststdev = NaN(nMonths,horizons,size(combo_specifications,1));
fcstskew = NaN(nMonths,horizons,size(combo_specifications,1));

%% loop across specifications 

tic;

for spec = 1:size(combo_specifications,1)
    
    % display spec
    fprintf('Starting loop iteration: spec = %d\n', spec);

    varnames = combo_specifications(spec,:); % variables you want to import 
    lagctryvar = ones(size(varnames)); % generate lagged variables and append them to the dataset 
    lagglobalvar = []; % lag for global var (if ygrowtrate = 1 then x2)
    lags = [lagctryvar lagglobalvar];
    
    %% Load the data
    
    dataloading_monthly; % load the data
    
    % get the estimation index
    idx_estimation = find(dateNumeric== StartEst); 
    idx_covid = find(dateNumeric== min(covid_dates)); 
    
    varnames_lag = strcat('l1', varnames); % 1-period lag to mimic BOE information set
    
    % dependent and explanatory variables
    depvarname  = varnames_lag(1,1); % Dependent variable name
    explvarname = varnames_lag(1,2:end); % Explanatory variable names

    if covid == 1
        explvarname = varnames_lag(1,2:end-1); % Explanatory variable names
        exovarname = dummyvarname; % Explanatory variable names % placeholder for exovarname = dummyvarname; 
    end

    %% Panel structure 
    
    % for now no panel reg
    if onlyuk == 1
    
        T = countryData.UK; %UK
        
        depvar = T.(depvarname{1}); % extract the dep var
        
        explvar = table2array(T(:, explvarname)); % extract the explanatory var

        if covid == 1
            exovar = T.(exovarname{1});
            exovar = [zeros(1, 1); exovar(1:end-1)]; % update to mimic info set of BOE
        end
   
    else
    
        % Extract panel data by stacking data from all countries
        countryNames = fieldnames(countryData);
        
        % Initialize empty arrays 
        depvar = [];
        explvar = [];
        
        % Loop over each country's table
        for i = 1:length(countryNames)
    
            T = countryData.(countryNames{i});
            
            % Extract dependent variable and explanatory variables for
            % country i
            dep_temp = T.(depvarname{1});
            expl_temp = table2array(T(:, explvarname));
            
            % Stack the data vertically
            depvar = [depvar; dep_temp];
            explvar = [explvar; expl_temp];
    
        end
    end

    %explvar_all_spec(:,:,spec) = explvar;
    %save(fullfile(outputFolder, 'explanatoryvar_infl.mat'),'explvar_all_spec', 'exovar'); % save explanatory variables and dummy covid 
    %meanVecExpl = mean(explvar, 1);   
    %stdVecExpl  = std(explvar, 0, 1);   
    %momentsVecExpl = [meanVecExpl;stdVecExpl];
    %momentsVecExpl = round(momentsVecExpl, 2);

    %% Regressions
    
    if model_selection == 1
    
    %%%%%%%%%%%%%%%%%%%%%
    %% OLS regression %%%
    %%%%%%%%%%%%%%%%%%%%%
    
    % pre-allocate the matrix 
    beta_hat = zeros(size(explvar,2)+1, size(depvar,1) - idx_estimation+1,horizons); % OLS coefficents 
    varcovbeta = zeros(size(explvar,2)+1,size(explvar,2)+1,  size(depvar,1) - idx_estimation+1,horizons); % var/cov of OLS coeff
    
        for endtime = idx_estimation:size(depvar,1)
        
            for h = 1:horizons % Loop over the horizon
                    
                % Generate the dependent variable
                dvar = depvar(h+1:endtime,:) ; 
                
                % generate the explanatory variable
                exvar = explvar(1:endtime-h,:);
            
                X = [ones(size(exvar,1), 1), exvar];
                y = dvar;
                
                % Run OLS regression
                beta_hatt = X \ y;
                beta_hat(:,endtime-idx_estimation+1,h) = beta_hatt;
        
                % Compute predicted values and residuals
                y_hat = X * beta_hatt;
                residuals = y - y_hat;
            
                % Calculate degrees of freedom
                n = size(X, 1);
                k = size(X, 2);
                
                % Estimate variance of the error terms 
                sigma2 = sum(residuals.^2) / (n - k);
            
                % get the variance covariance matrix of the coefficent 
                varcovbeta(:,:,endtime-idx_estimation+1,h) = sigma2 * inv(X'*X);
            
            end
        
        end
        
        % get the forecast mean and the forecast variance 
        for endtime = idx_estimation:size(depvar,1)

            for h = 1:horizons % Loop over the horizon
        
                fcstmean(endtime-idx_estimation+1,h,spec) =  [1, explvar(endtime,:)] * beta_hat(:,endtime-idx_estimation+1,h); % mean fcst
                fcststdev(endtime-idx_estimation+1,h,spec) = sqrt([1, explvar(endtime,:)] * varcovbeta(:,:,endtime-idx_estimation+1,h) * [1, explvar(endtime,:)]'); % sd fcst
                % skewness is empty
            end
           
        end     
      
    elseif model_selection == 2 || model_selection == 3 || model_selection == 4
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Quantile regression %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for endtime = idx_estimation:size(depvar,1)

                % display spec
                fprintf('loop over time: %d out of %d\n', endtime - idx_estimation +1, size(depvar,1) - idx_estimation + 1);
                                
                % Generate the dependent variable
                dvar = depvar(1:endtime,:) ;
           
                % generate the explanatory variables
                X = explvar(1:endtime,:);
                
                % store the dep variable in a comformable 3d matrix
                Y_LP = NaN (size(X,1),1,horizons);

                for hor=1:horizons

                    Y_LP(1:size(X,1)-hor,1,hor) = dvar(1+hor:end);
                 
                end

                if covid == 0 % NO COVID

                    if endtime ~=size(depvar,1) % if we are not in the last period -> no bootstrap
                        
                       bst = 0;
                    
                    else % if we are in the last period -> bootstrap
                        
                       bst = 1;

                    end

                    [model.bQR,model.bQRbst,~,~] = qfe_qr_local_projection_SL_final(Y_LP,[ones(size(X,1),1), X],quantilelevels,(1:horizons)',bst,bstOptions, 0, 0); % qreg
                    coeffqr_SL(1:end-1,:,:,endtime-idx_estimation+1,spec) = model.bQR; % store the coeff
                    bootstrapqrg_SL(1:end-1,:,:,:,spec)                  = model.bQRbst; % store the bootstrap   

                    % loop over horizon
                    for h = 1:horizons % Loop over the horizon
                        predicted_quantiles_SL(endtime-idx_estimation+1,:, h, spec) = [1, X(endtime,:), 0] * coeffqr_SL(:,:,h, endtime-idx_estimation+1,spec);
                    end

                else  % COVID

                    if endtime <= idx_covid % pre-covid -> no boostrap

                        covid_qreg = 0; % no covid
                        exo = 0; % empty exo
                        bst = 0; % no boots
                        place_holder = size(X,2) + 1; % for the constant 

                    else  % post-covid
    
                        covid_qreg = 1; % covid dummy
                        exo = exovar(1:endtime,:); % dummy covid    
                        place_holder = size(X,2) + 2; % for covid and the constant

                        if endtime < size(depvar,1) % post-covid and not in the final period -> no boostrap

                            bst = 0;

                        elseif endtime == size(depvar,1) % post-covid and final period -> boostrap

                            bst = 1;

                        end

                     end
                    

                    [model.bQR,model.bQRbst,~,~] = qfe_qr_local_projection_SL_final(Y_LP,[ones(size(X,1),1), X],quantilelevels,(1:horizons)',bst,bstOptions,covid_qreg, exo);
                    coeffqr_SL(1:place_holder,:,:,endtime-idx_estimation+1,spec) = model.bQR; % store the coeff
                    bootstrapqrg_SL(1:place_holder,:,:,:,spec)                  = model.bQRbst;   % store the bootstrap  
                    
                    % loop over horizon
                    for h = 1:horizons % Loop over the horizon
                        predicted_quantiles_SL(endtime-idx_estimation+1,:, h, spec) = [1, X(endtime,:), exo(end)] * coeffqr_SL(:,:,h, endtime-idx_estimation+1,spec);
                    end

                end

        end   

        % sort the quantiles
        predicted_quantiles_SL = sort(predicted_quantiles_SL, 2);     

    end

    if model_selection == 2
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %% Skew T Distribution %%
        %%%%%%%%%%%%%%%%%%%%%%%%%

        % get the forecast mean and the forecast variance

        for endtime = idx_estimation:size(depvar,1)

        fprintf('loop over time: %d out of %d\n', endtime - idx_estimation +1, size(depvar,1) - idx_estimation + 1);

            for h = 1:horizons % Loop over the horizon

                [lc_skt(endtime-idx_estimation+1,h,spec), sc_skt(endtime-idx_estimation+1,h,spec), sh_skt(endtime-idx_estimation+1,h,spec), df_skt(endtime-idx_estimation+1,h,spec), ~, ~, ~, ~] = QuantilesInterpolation(predicted_quantiles_SL(endtime-idx_estimation+1,:,h,spec), quantilelevels);
                [fcstmean(endtime-idx_estimation+1, h, spec), fcststdev(endtime-idx_estimation+1, h, spec), fcstskew(endtime-idx_estimation+1, h, spec), ~] = moments_skewt_distr_sm_updated(lc_skt(endtime-idx_estimation+1,h,spec), sc_skt(endtime-idx_estimation+1,h,spec), sh_skt(endtime-idx_estimation+1,h,spec),df_skt(endtime-idx_estimation+1,h,spec));
           
            end

        end

    elseif model_selection == 3
    
        %%%%%%%%%%%%%%%%%%%%%%
        %% Semi- parametric %%
        %%%%%%%%%%%%%%%%%%%%%%

        % get the empirical moments 

        for h = 1:horizons % Loop over the horizon
    
            semi_param_distr(:,:,h,spec) = QR_sm(predicted_quantiles_SL(:,:,h,spec), quantilelevels); % get the semi-param distr
            fcstmean(:,h,spec) = mean(semi_param_distr(:,:,h,spec),2); % mean 
            fcststdev(:,h,spec) = std(semi_param_distr(:,:,h,spec), 0, 2); % sd 
            fcstskew(:,h,spec)  = skewness(semi_param_distr(:,:,h,spec), 0, 2); % skew

        end
        
        % get the empirical cdf

          for h = 1:horizons % Loop over the horizon

            for endtime = idx_estimation:size(depvar,1)

                    [f,x] = ecdf(semi_param_distr(endtime-idx_estimation+1,:,h,spec));
                    block     = nan(20001,2);                 % fixed-size holder for handling duplicates
                    block(1:numel(x),:) = [x,f];              % copy what we got
                    empirical_cdf_semi_param(:,:,endtime-idx_estimation+1, h,spec) = block; % cdf

            end      
        
          end

    elseif model_selection == 4

        %%%%%%%%%%%%%%%%%%%%%%
        %% Two-piece Normal %%
        %%%%%%%%%%%%%%%%%%%%%%

    %load(fullfile(outputFolder, 'predicted_quantiles_inflation.mat'),'predicted_quantiles_SL');  % across all specifications (time x quantile x hor x spec)
    %for spec = 1:size(combo_specifications,1)
        % get the forecast mean and the forecast variance and fcst skew
        for endtime = idx_estimation:size(depvar,1)

            fprintf('loop over time: %d out of %d\n', endtime - idx_estimation +1, size(depvar,1) - idx_estimation + 1);
     
            for h = 1:horizons % Loop over the horizon
    
                param_tpn(endtime-idx_estimation+1, :, h, spec) = fit_tpn_to_quantiles_SM(predicted_quantiles_SL(endtime-idx_estimation+1,:,h,spec), quantilelevels);
                
                fcstmean(endtime-idx_estimation+1,h,spec) = two_part_normal_mean(param_tpn(endtime-idx_estimation+1, :, h, spec)) ;
                fcststdev(endtime-idx_estimation+1,h,spec) = sqrt(two_part_normal_variance(param_tpn(endtime-idx_estimation+1, :, h, spec))) ;
                fcstskew(endtime-idx_estimation+1,h,spec) = two_part_normal_skewness(param_tpn(endtime-idx_estimation+1, :, h, spec)) ;
                
            end
               
        end
   %end

    end

end

toc;

% save coeff - bootstrapped coeff, predicted quantiles and moments
%save(fullfile(outputFolder, 'qreg_results_inflation.mat'),'coeffqr_SL');  % across all specifications (variables x quantile x hor x time x spec) 
%save(fullfile(outputFolder, 'predicted_quantiles_inflation.mat'),'predicted_quantiles_SL');  % across all specifications (time x quantile x hor x spec)
%save(fullfile(outputFolder,'bootstrap_results_inflation.mat'),'bootstrapqrg_SL'); % across all the specifications but only for the last period  (variables x quantile x hor x bootnumb x spec)  (h = 13)

% save empirical cdf for the best model
% save(fullfile(outputFolder,'semi_param', 'moments_all_spec.mat'),'fcstmean', 'fcstskew', 'fcststdev'); % across all the specifications mean sd and skew time x hor x spec
%empirical_cdf_semi_param_best = empirical_cdf_semi_param(:,:,:,:,idx_min);
%save(fullfile(outputFolder, 'semi_param', 'empirical_cdf_best_model_semi_param.mat'),'empirical_cdf_semi_param_best', '-v7.3');  % sample x 2 (support of the cdf and values of the cdf) x time x horizons (only for the best model) [40 sec to save]

% save skew t parameters and moments
%save(fullfile(outputFolder,'sktparam', 'skewtparam_all_spec.mat'),'lc_skt', 'sc_skt', 'sh_skt', 'df_skt'); % across all the specifications mean sd and skew time x hor x spec
%save(fullfile(outputFolder,'sktparam', 'moments_all_spec.mat'),'fcstmean', 'fcstskew', 'fcststdev'); % across all the specifications mean sd and skew time x hor x spec
%lc_skt_best = lc_skt(:,:,idx_min); sc_skt_best = sc_skt(:,:,idx_min); sh_skt_best = sh_skt(:,:,idx_min); df_skt_best = df_skt(:,:,idx_min);
%save(fullfile(outputFolder,'sktparam', 'skewtparam_best_spec.mat'),'lc_skt_best', 'sc_skt_best', 'sh_skt_best', 'df_skt_best'); % across all the specifications mean sd and skew time x hor x spec

% two piece normal 
%save(fullfile(outputFolder,'two_piece_normal', 'moments_all_spec.mat'),'fcstmean', 'fcstskew', 'fcststdev'); % across all the specifications mean sd and skew time x hor x spec
%save(fullfile(outputFolder,'two_piece_normal', 'tpn_all_spec.mat'),'param_tpn'); % across all the specifications mode sd1 and sd2
%param_tpn_best = param_tpn(:,:,:,idx_min);
%save(fullfile(outputFolder,'two_piece_normal', 'tpn_best_spec.mat'),'param_tpn_best'); % across all the specifications mode sd1 and sd2


%%%%%%%%%%
%% CRPS %%
%%%%%%%%%%

%% Import acutal data (only for UK)
actual_var_long = countryData.(ctrynames{1}).(dep_var_name{1}); 

% save actual inflation 
% save(fullfile(outputFolder, 'actual_infl.mat'),'actual_var_long'); % save
% actual inflation

%covid_dummy_infl =  countryData.(ctrynames{1}).covid; 
%save(fullfile(outputFolder, 'covid_dummy_infl.mat'),'covid_dummy_infl'); % save explanatory variables and dummy covid 

% Filter table for the estimation period 
actual_var_long_filtered = actual_var_long(idx_estimation:end, :); % note: idx_estimation come from the loop but it is always the same as it is linked to StartEst in the control panel

% Number of observations
n = length(actual_var_long_filtered);

% Initialize matrix to hold actual data
interm = NaN(n, n);

% Fill the matrix with the required observations
for i = 1:n
    interm(1:(n-i+1), i) = actual_var_long_filtered(i:end);
end

actualvar = interm(1:horizons,:)'; % final product of actual data (note that n has to be greater than horizons)

% per-allocate vector of names
spec_names = strings(size(combo_specifications,1),1);

% construct names of the specifications
for spec = 1:size(combo_specifications,1)
    spec_names(spec) = strjoin(combo_specifications(spec,:), ' ');
end

%% Run CRPS
    
if model_selection == 3

    % pre-allocate the matrix (estimated time per specification - 
    crps_results = zeros(size(semi_param_distr,1),horizons,size(semi_param_distr,4));
    
    for spec = 1:size(semi_param_distr,4)
    
        for h = 1:horizons
    
            for time = 1: size(semi_param_distr,1)
    
            crps_results(time,h,spec)  = crps(semi_param_distr(time,:,h,spec),actualvar(time,h),2 ); % run the crps
    
            end
    
        end
    end
    
    % save matlab objects
    % save(fullfile(outputFolder, 'crps', 'results_inflation.mat'), 'crps_results');
  
    % we evaluate the crps in the first 2 year 
    crps_results = crps_results(:, 1:25, :); 

    % get the average crps for each specification
    avg_crps_x_hor = squeeze(mean(crps_results, 1,  'omitnan')); % avg over time for a given hor and spec 
    avg_crps = squeeze(mean(crps_results, [1 2], 'omitnan')); % avg over time and over hor for a given spec
    
    % pre allocate by hor results 
    idx_min_hor = zeros(numel(horizons_crps),1);
    
    for i = 1:numel(horizons_crps) 
    
        hor = horizons_crps(i);
    
        select_row = avg_crps_x_hor(hor, :); % row of interest
        
        [~, idx_min_hor(i)] = min(select_row); % best x hor
    
    end
    
    % best specification across all hor
    [~, idx_min] = min(avg_crps);
    fprintf('The best model is:\n  %s\n(with CRPS = %.4f)\n\n', spec_names(idx_min), avg_crps(idx_min));
    
    %  SAVE TO EXCEL: sheet 1 is avg over time and over hor for a given spec 
    %                 sheet 2 is best spec per horizon 
    
    % sheet 1 
    T1 = table(spec_names, avg_crps, ...
              'VariableNames', {'Specification','CRPS'});
    writetable(T1, fullfile(outputFolder, 'crps', sprintf('crps_results_%s_inflation.xlsx', modellist{model_selection})), 'Sheet', 'best_spec', 'WriteMode','overwrite')
    
    % build a table listing each horizon and its best spec name
    bestSpecNames_hor = spec_names(idx_min_hor);
    T2 = table(horizons_crps(:), bestSpecNames_hor, ...
               'VariableNames', {'Horizon','BestSpecification'});
    % write to a second sheet called “best_by_horizon”
    writetable(T2, fullfile(outputFolder, 'crps', sprintf('crps_results_%s_inflation.xlsx', modellist{model_selection})), 'Sheet', 'best_spec_by_hor', 'WriteMode','append');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE THE MOMENTS INTO AN EXCEL FILE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% condition that tells which specification to save
if use_best_spec == 1 
    
    spec_to_use = idx_min; %idx
    
elseif use_best_spec == 0
    
    spec_to_use = specplot;
        
end

% Generate quarterly dates
startDate_est = datetime(StartEst, 'ConvertFrom', 'datenum');  
endDate = datetime(endT, 'ConvertFrom', 'datenum');
    
% Create a vector of dates
months_date = startDate_est:calmonths(1):endDate;
   
if use_best_spec == 1 
        
    % Define Excel filename once using the model name
    filename = fullfile(outputFolder, sprintf('%s_best_spec.xlsx', modellist{model_selection}));
else
    % Define Excel filename once using the model name
    filename = fullfile(outputFolder, sprintf('%s_spec_%d.xlsx', modellist{model_selection}, specplot));
end

% Loop over each variable in the list
for i = 1:length(momentlist)
    
%    % Get the current forecast data
    currentVar = eval(momentlist{i});
    currentVar_final = currentVar(:,:,spec_to_use);
        
%    % Create the table with quarterly dates as the first column
    T = [table(months_date', 'VariableNames', {'Dates'}), array2table(currentVar_final)];
        
    % Rename the forecast columns from the second column onward as h_0, h_1, ..., h_12
    forecastNames = strcat("h_", string(0:horizons-1));
    T.Properties.VariableNames(2:end) = cellstr(forecastNames);
        
    % Export the table to a specific sheet in the Excel file
    writetable(T, filename, 'Sheet', momentlist{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IMPULSE RESPONSE FUNCTION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load(fullfile(outputFolder, 'qreg_results_inflation.mat'), 'coeffqr_SL');
%coeffqreg_SL = permute(coeffqr_SL, [2 1 4 3 5]); % quant x var x time x hor x spec
%load(fullfile(outputFolder, 'bootstrap_results_inflation.mat'), 'bootstrapqrg_SL');
%load(fullfile(outputFolder, 'analysis_2000_2019', 'qreg_results_inflation.mat'), 'coeffqr_SL');
%coeffqreg_SL = permute(coeffqr_SL, [2 1 4 3 5]); % quant x var x time x hor x spec
%load(fullfile(outputFolder , 'analysis_2000_2019', 'bootstrap_results_inflation.mat'), 'bootstrapqrg_SL');
% coeffqreg_SL(:,:,:,14:end,:) = [];
%bootstrapqrg_SL(:,:,14:end,:,:) = [];

if model_selection ~= 1

    % condition that tells which specification to plot
    if use_best_spec == 1 
    
        spec_to_use = idx_min; %idx
        varnames_to_plot = combo_specifications(idx_min,2:end); % names
    
    elseif use_best_spec == 0
    
        spec_to_use = specplot;
        varnames_to_plot = combo_specifications(specplot,2:end);
        
    end
    
    % get the index of the quantiles you want to plot from quantilelevls plot 
    idx_qt = arrayfun(@(q) find(abs(quantilelevels - q) < 1e-8, 1, 'first'), quantilesirf);
    
    % choose one color per quantile and make the legend labelicer
    colors = lines(numel(idx_qt)); 
    labels = arrayfun(@(q) sprintf('%.0f^{th}',q*100), quantilesirf, 'UniformOutput', false);
    
    % chart
    figure('Position',[100 100 1200 800], 'Color','w');
    
    t = tiledlayout(3,2, 'TileSpacing','Compact', 'Padding','Compact'); % always 5 expl vars so 3x2 
    
    % Preallocate for legend 
    h = gobjects(numel(idx_qt),1);
    
    for ii = 1:numel(varnames_to_plot) % loop over variables
        
        ax = nexttile;       
        hold(ax,'on');
        
        for dd = 1:numel(idx_qt) % loop over quantiles
    
            q_idx    = idx_qt(dd);
            % coeff_irf = squeeze( coeffqreg(q_idx, 1+ii, end, :, spec_to_use) );
            % coeff_irf = squeeze( coeffqreg_SL(q_idx, 1+ii, end, :, spec_to_use) );
   
            % bootstrap draws: horizons×nBoot
            % bs_draws = squeeze( bootstrapqrg(1+ii, q_idx, :, :, spec_to_use) );
            bs_draws = squeeze( bootstrapqrg_SL(1+ii, q_idx, :, :, spec_to_use) );
            coeff_irf = prctile(bs_draws, 50, 2);   % horizons×1
          
            % compute 25th and 75th percentiles at each horizon
            lb = prctile(bs_draws, 16, 2);   % horizons×1
            ub = prctile(bs_draws, 84, 2);   % horizons×1
            %lb = squeeze(lb_qr(1+ii,q_idx,:));
            %ub = squeeze(ub_qr(1+ii,q_idx,:));
    
            if ii == 1 % store info for the legend from the first chart
                h(dd) = plot(ax, 1:horizons, coeff_irf, 'LineWidth',1.5, 'Color', colors(dd,:) );
                if plot_bands_irf == 1
                    fill(ax, [1:horizons, fliplr(1:horizons)], [ub', fliplr(lb')], ...
                    colors(dd,:), 'FaceAlpha', 0.1, 'EdgeColor','none');
                end
            else
                plot(ax, 1:horizons, coeff_irf, 'LineWidth',1.5, 'Color', colors(dd,:) );
                if plot_bands_irf == 1
                    fill(ax, [1:horizons, fliplr(1:horizons)], [ub', fliplr(lb')], ...
                    colors(dd,:), 'FaceAlpha', 0.1, 'EdgeColor','none');
                end 
            end
        end
        
        yline(ax, 0, 'k--', 'LineWidth',1); % zero line
        grid(ax,'on');
        xlim(ax,[1 horizons]);
        xlabel(ax,'Quarters','FontSize',10);
        ylabel(ax,'','FontSize',10);
        title(ax, varnames_to_plot{ii}, 'FontSize',12,'FontWeight','bold');
        ax.FontSize = 10;
        hold(ax,'off');
    
    end
    
    % Create a single, shared legend from the handles
    lg = legend(h, labels,'Orientation', 'horizontal', 'Location', 'northoutside', 'Box', 'off', 'FontSize',  12);
    lg.Layout.Tile = 'north';
    
    % Overall title for the layout
    title(t, sprintf('Model: %s \nSpecification: %s', modellist{model_selection}, spec_names{spec_to_use}), 'FontSize',14);
    
    % save the title
    %fileName = sprintf('irfs_end_period_%s_model_spec_%s_bands_h13_median', modellist{model_selection}, spec_names{spec_to_use});
    %saveas(gcf, fullfile(outputFolder, 'irfs', [fileName, '.png']));

end

%% Other potential options (to add in the future)
% xGFC = 0;
% Standardise y variable? (1 if yes, 0 if no) [default = 0]
% standardy   = 0;
% Cumulate growth of y in local projection? (1 if yes, 0 if no)
% cumulateY   = 0;
% Average annual growth of y variable? (1 if yes, 0 if no) [default = 1]
% opt_avg     = 0;
% Which to standardise? (1 if yes, 0 if no)
% varDomstd    = [zeros(1,length(ctryList)),0,1,0];    
% Select if there are foreign variables
%forwts = 1;     % 1 if yes, 0 if no
% orthog_domvars = 1; %1 if orthogonalise domestic variables wrt foreign variables, 0 otherwise (NB switches to zero if no foreign variables)
% Note: potentially you can change the weigths 

% OTHERS

% bootstrapped OLS

%{

horizon = 1:horizons;
p25 = prctile(bootols,5,3);   % 6×13
p75 = prctile(bootols,95,3);   % 6×13
varnames_to_plot = combo_specifications(1,2:end); % names

figure;
for i = 1:5
    subplot(3,2,i);
    hold on;

    % shaded band
    xpoly = [horizon, fliplr(horizon)];
    ypoly = [p25(i+1,:), fliplr(p75(i+1,:))];
    fill(xpoly, ypoly, [0.8 0.8 0.8], ...
         'EdgeColor','none','FaceAlpha',0.4);

    % OLS line
    plot(horizon, coeffols(i+1,:), 'k-', 'LineWidth', 2);

    % optional percentile lines
    plot(horizon, p25(i+1,:), 'r--', 'LineWidth', 1);
    plot(horizon, p75(i+1,:), 'r--', 'LineWidth', 1);

    title(varnames_to_plot{i});
    xlabel('Horizon');
    ylabel('Coefficient');
    grid on;
    hold off;
end

% if you want an overall title
sgtitle('OLS Coefficients with 5–95% Bootstrap Bands');

%}

%% REPLICATION FIGURE 1 pag 5 - LOPEZ SALIDO AND LORIA (2024)

% condition that tells which specification to plot
if use_best_spec == 1 

    spec_to_use = idx_min; %idx
    % varnames_to_plot = combo_specifications(idx_min,2:end); % names
    varnames_to_plot = { sprintf('INFLATION\nPERSISTENCE'),  sprintf('INFLATION\nEXPECTATIONS'), ...
                         sprintf('ECONOMIC\nTIGTHNESS'), sprintf('EXTERNAL\nCONDITIONS'), ...
                         sprintf('FINANCIAL\nCONDITIONS')};  % uncomment it if you run other variables  (, ...sprintf('DUMMY\nCOVID'))

elseif use_best_spec == 0

    spec_to_use = specplot;
    varnames_to_plot = combo_specifications(specplot,2:end);
    
end

load(fullfile(outputFolder, 'qreg_results_inflation.mat'), 'coeffqr_SL'); % load coeff - for CPI is _CPI - _CPI_no_ener
coeffqreg_SL = permute(coeffqr_SL, [2 1 4 3 5]); % quant x var x time x hor x spec
load(fullfile(outputFolder, 'bootstrap_results_inflation.mat'), 'bootstrapqrg_SL'); % load bootstrap
% load(fullfile(outputFolder, 'analysis_2000_2019', 'bootstrap_results_inflation.mat'), 'bootstrapqrg_SL'); % load bootstrap


quantilesirf = [0.05 0.1 0.25 0.5 0.75 0.90 0.95]; % quantiles in the chart

% get the index of the quantiles you want to plot from quantilelevls plot 
idx_qt = arrayfun(@(q) find(abs(quantilelevels - q) < 1e-8, 1, 'first'), quantilesirf);

% select the horizon of interest 
hor_of_interest = [2, 13, 25]; % 1m-2-3y ahead

% choose one color per quantile and make the legend labelicer
colors = lines(numel(idx_qt)); 
labels = arrayfun(@(q) sprintf('%.0f^{th}',q*100), quantilesirf, 'UniformOutput', false);

% compute the standard deviation of the series
% median_qr   = squeeze(median(bootstrapqrg_SL(:,:,:,:,spec_to_use),4)); % this is the median over the bootstraps so the size is # of variables x # of quantiles x horizon

% Calculate the standard deviations for each variable, quantile and horizon over bootstraps
std_qreg = squeeze( std( bootstrapqrg_SL(:,:,:,:,spec_to_use) ,0,4) );

% chose the order of the variables (to compare with Lopez-Salido)
var_order = [1, 2, 3, 5, 4];  % infl persistence - infl exp - econ slack - financial market - external sector (add 6 to check dummy covid)

%figure('Position',[100 100 900 1200],'Color','w'); % size
fig = figure('Units','centimeters','Position',[1 1 22 30],'Color','w');

t = tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');

for k = 1:numel(var_order) % var 

    ii = var_order(k); % access the series you want to plot - in the order you want to plot it
    
    for iH = 1:numel(hor_of_interest) % hor

        hor = hor_of_interest(iH); % get the horizon
        ax = nexttile((k-1)*numel(hor_of_interest) + iH); % logic: fix a variable and then plot all the horizons. Then move to the next variable and so on..
    
        hold(ax,'on');
        yline(ax,0,'k'); % zero line
        
        % extract point estimate 
        coeffs = nan(1,numel(idx_qt)); % pre-allocate the vectors
        sds = nan(1,numel(idx_qt)); % pre-allocate the vectors
    
        for dd = 1:numel(idx_qt)
    
            q_idx        = idx_qt(dd); % get the index of the quantile
            coeffs(dd)   = squeeze( coeffqreg_SL(q_idx,1+ii,end,hor,spec_to_use) ); % qtl x var x hor x spec (last period is used)
            % coeffs(dd)     = median_qr(1+ii,q_idx,hor);
            sds(dd)       = std_qreg(1+ii,q_idx,hor);
            
            % Plot marker‐only errorbar 
            errbar = errorbar(ax, dd, coeffs(dd), sds(dd), 's', ...
                         'LineWidth',1.5, 'CapSize',0); 
            errbar.MarkerFaceColor = colors(dd,:);
            errbar.MarkerEdgeColor = 'none';
            
        end
        
    % Axis formatting
    ax.XLim       = [0.5, (numel(idx_qt)+0.5)]; 
    ax.XTick      = 1:numel(idx_qt);
    ax.XTickLabel = labels;
    
    % y axis
    if k == 1  % inflation persistence
        ylim([-2,  1]);
    elseif k == 2  % inflation expectations
        ylim([-1,  4]);
    elseif k == 3  % 1y delta v/u
        ylim([-15,  15]);
    elseif k == 4  % bond spread
        ylim([-1.5e-2,  1.5e-2]);
        ax.YAxis.Exponent = -2;   % force ×10⁻² on the y-axis
    elseif k == 5  % oil inflation
        ylim([-0.20, 0.4]);
        ax.YAxis.Exponent = -2;   % force ×10⁻² on the y-axis
    end
    ax.YAxis.TickLabelFormat = '%.1f';

    if iH==1  % title hor only for top charts
      ylabel(ax, varnames_to_plot{ii}, 'Interpreter','none');
    end
    if k==1 % title var only for the colomn on the left
    
        if hor == 2
            
            title(ax, sprintf('%g-MONTH', 1), 'Units', 'normalized', 'Position', [0.5, 1.10, 0]);

        else
           
            value = (hor - 1) / 12;  
            title(ax, sprintf('%g-YEAR', value), 'Units', 'normalized', 'Position', [0.5, 1.10, 0]);
        
        end

    end

    ax.FontSize        = 10;   % ticks & axis labels
    ax.Title.FontSize  = 12;   % subplot titles

    box(ax,'on');
    hold(ax,'off');
    
    end 
end

% Overall title
%title(t, sprintf('Model: %s  \nSpec: %s', modellist{model_selection}, spec_names{spec_to_use}), ...
%       'FontSize',10);

set(fig, 'PaperUnits','centimeters','PaperPosition',[0 0 22 30]);

% Save as high-res PNG
%print(fig, fullfile(dropboxFolder, 'econ_interpretation_charts', ...
%    'econ_interpr_inflation_1m_1y_2y_no_dcovid.png'), '-dpng', '-r300');

% save the title
% saveas(gcf, fullfile(dropboxFolder, 'econ_interpretation_charts', 'econ_interpr_inflation_1m_1y_2y.png'));

% plot evolution of coefficents over time 
%hor_of_interest = [2, 13, 25];
%var_of_interest = 4;

% Dates (quarterly from 2004q1 to 2025q1)
%dates = datetime(2004,1,1):calmonths(1):datetime(2025,1,1);

% Plot for each horizon
%for i = 1:numel(hor_of_interest)
%    h = hor_of_interest(i);
    
%    figure; hold on;
%    plot(dates, squeeze(coeffqreg_SL(2,  var_of_interest, :, h, spec_to_use)), 'LineWidth', 1.5);
%    plot(dates, squeeze(coeffqreg_SL(10, var_of_interest, :, h, spec_to_use)), 'LineWidth', 1.5);
%    plot(dates, squeeze(coeffqreg_SL(18, var_of_interest, :, h, spec_to_use)), 'LineWidth', 1.5);
    
    % Legend above chart
%    legend('10th','50th','90th', ...
%        'Location','northoutside','Orientation','horizontal');
%    xlabel('Date'); ylabel('Coefficient');
%    legend('10th','50th','90th');
%    grid on;
%end

%% REPLICATION FIGURE 2 pag 6 - LOPEZ SALIDO AND LORIA (2019)

% condition that tells which specification to plot
if use_best_spec == 1

    spec_to_use = idx_min; %idx

elseif use_best_spec == 0

    spec_to_use = specplot;

end

% load the predicted quantiles

% out-of-sample forecasts
load(fullfile(outputFolder, 'predicted_quantiles_inflation.mat'),'predicted_quantiles_SL');  % recursive - across all specifications (time x quantile x hor x spec) for CPI is _CPI

% in-sample forecasts 
%load(fullfile(outputFolder, 'qreg_results_inflation.mat'), 'coeffqr_SL'); % qreg 
%coeffqreg_SL = permute(coeffqr_SL, [2 1 4 3 5]); % reshape qreg -> quant x var x time x hor x spec 
%coeffqreg_last = squeeze(coeffqreg_SL(:,:,end,:)); % consider only the last period 

%load(fullfile(outputFolder, 'explanatoryvar_infl.mat'),'explvar_all_spec', 'exovar'); % load explanatory variables acrosss all specifications and dummy covid 
%predictors = [ones(size(exovar)), explvar_all_spec, zeros(size(exovar))]; % predictors = [ones(size(exovar)), explvar_all_spec, exovar]; % matrix of predictors 
%predictors = predictors(168:end,:); % only for the period of interest (2003m12 onward --> first forecast is 2004m1 1 period ahead)

% get predicted quantiles
%for hh = 1:horizons
%    predicted_quantiles_SL(:,:,hh) = predictors * coeffqreg_last(:,:,hh)' ;
%end
        
%predicted_quantiles_SL = sort(predicted_quantiles_SL, 2);  % sort the quantiles 

% cross check last value obtained with in sample is same as oos

load(fullfile(outputFolder, 'actual_infl.mat'),'actual_var_long'); % load actual infl for the full estimation period 1990m1-2025m1 - notice that if we change the est period we need to change this variable - for CPI is _CPI
% load(fullfile(outputFolder, 'covid_dummy_infl.mat'),'covid_dummy_infl'); % load covid dummy
idx_estimation = 168; % 2004m1 is 168

% convert to datetime to construct the date variable for the chart
StartEstDT = datetime(StartEst, 'ConvertFrom','datenum'); % first period at which fcst are estimated
endTDT =  datetime(endT, 'ConvertFrom','datenum'); % last period at which the forecast are performed

quantilesplot = [0.05 0.10 0.25 0.50 0.75 0.90 0.95];
idx_qt = arrayfun(@(q) find(abs(quantilelevels - q) < 1e-8, 1, 'first'), quantilesplot);
bands = {'5^{th}–95^{th}','10^{th}–90^{th}','25^{th}–75^{th}'};

% get actual inflation (168 is m1-2004)
actual_var_long_filtered = actual_var_long(idx_estimation:end, :); 
% covid_dummy_filtered = covid_dummy_infl(idx_estimation:end, :); 

h_list = [2, 4, 13]; % horizon that we want to plot

for ii = 1:numel(h_list)

    h_plot = h_list(ii);      % horizon of interest   
    
    inf_to_plot_evol = actual_var_long_filtered(h_plot:end); % shift inflation forward to match the forecast
    
    % prepare the shaded area 
    % dummy_covid_to_plot = covid_dummy_filtered(1:end-h_plot); % shift to match the forecast

    months = (StartEstDT + calmonths(h_plot-1)) : calmonths(1) : endTDT;       % get the quaretrs 
    
    % get the quantile of interest 
    Q05 = squeeze(predicted_quantiles_SL(1:size(months,2), idx_qt(1), h_plot, spec_to_use));
    Q10 = squeeze(predicted_quantiles_SL(1:size(months,2), idx_qt(2), h_plot, spec_to_use));
    Q25 = squeeze(predicted_quantiles_SL(1:size(months,2), idx_qt(3), h_plot, spec_to_use));
    Q50 = squeeze(predicted_quantiles_SL(1:size(months,2), idx_qt(4), h_plot, spec_to_use));  % median
    Q75 = squeeze(predicted_quantiles_SL(1:size(months,2), idx_qt(5), h_plot, spec_to_use));
    Q90 = squeeze(predicted_quantiles_SL(1:size(months,2), idx_qt(6), h_plot, spec_to_use));
    Q95 = squeeze(predicted_quantiles_SL(1:size(months,2), idx_qt(7), h_plot, spec_to_use));
    
    figure('Position',[100 100 1200 800],'Color','w'); % size and white background
    
    tl = tiledlayout(1,1,'TileSpacing','compact','Padding','compact'); 
    ax = nexttile;

    hold on; 

    % 5–95 band (lightest blue)
    fill([months, fliplr(months)], [Q05', fliplr(Q95')], ...
         [0.85 0.9 1], 'EdgeColor','none','FaceAlpha', 1, 'DisplayName',bands{1});

    % 10–90 band (medium blue)
    fill([months, fliplr(months)], [Q10', fliplr(Q90')], ...
         [0.65 0.8 1], 'EdgeColor','none','FaceAlpha', 1, 'DisplayName',bands{2});

    % 25–75 band (darkest blue)
    fill([months, fliplr(months)], [Q25', fliplr(Q75')], ...
         [0.4 0.6 1], 'EdgeColor','none','FaceAlpha', 1, 'DisplayName',bands{3});
   
    % boundary lines for 5th and 95th percentiles = with dots
    plot(months, Q05, ':', 'Color',[0 0 0.7], 'LineWidth',1,'HandleVisibility', 'off');
    plot(months, Q95, ':', 'Color',[0 0 0.7], 'LineWidth',1,'HandleVisibility', 'off');

    yline(0, 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off'); % zer line

    yl = ylim;  % grab them for our patch

    % plot the shaded area
    %idx1 = find(dummy_covid_to_plot==1); % find the first occurrence of covid

    %if ~isempty(idx1)

    %    startIdx = idx1(1); % forecast done at 2020q2
    %    endIdx   = idx1(end); % forecast done at 2022q2
    %    ms = months(startIdx);
    %    me = months(endIdx);
    %    yl = ylim;  % get current y‑limits
        
        % draw exactly one grey patch
    %    patch([ms me me ms], ...
    %          [yl(1) yl(1) yl(2) yl(2)], ...
    %          [0.6 0.6 0.6], 'EdgeColor','none', ...
    %          'FaceAlpha',0.3, 'HandleVisibility','off');
    %    uistack(findobj(gca,'Type','patch'),'top'); % change top or bottom
    %end

    % median line
    plot(months, Q50, '--', 'Color', [0 0 0.7], 'LineWidth', 1, 'DisplayName','50^{th}');

    plot(months, inf_to_plot_evol, 'k','DisplayName','Outturn','LineWidth',1.25); % actual var
    
    % ylim([0, inf]);
    grid on;

    ticks = months(1):calyears(2):months(end); % every 2 years
    ax.XTick = ticks;
    ax.XAxis.TickLabelFormat = 'yyyy';   % show only the year
    
    % legend above the axes
    lgd = legend(ax,'show','Location','northoutside','Orientation','horizontal');
    lgd.Box = 'off';
    
    set(gca,'FontSize',18);
    set(lgd,'FontSize',18);

    hold off;

    %filename = sprintf('predict_distr_infl_h%d_no_dcovid.png', h_plot);
    %filepath = fullfile(dropboxFolder, 'predictive_densities', filename);
    %saveas(gcf, filepath);

end

%% HISTORICAL DECOMPOSITION CHART 

% note that here we assume all the mat element refers to the best (and
% only) spec --> spec dimension is 1 

% load the predicted quantiles
load(fullfile(outputFolder, 'predicted_quantiles_inflation.mat'),'predicted_quantiles_SL');  % recursive - across all specifications (time x quantile x hor x spec) for CPI is _CPI

load(fullfile(outputFolder, 'actual_infl.mat'),'actual_var_long'); % load actual infl for the full estimation period 1990m1-2025m1 - notice that if we change the est period we need to change this variable - for CPI is _CPI
load(fullfile(outputFolder, 'covid_dummy_infl.mat'),'covid_dummy_infl'); % load covid dummy
idx_estimation = 168; % 2004m1 is 168
spec_to_use = 1;

varnames_plot = [{'cons'}, combo_specifications(2:end), {'dcovid'}]; % varnames
varnames_plot = cellfun(@(x) strrep(x,'_',' '), varnames_plot, 'UniformOutput', false);

% in-sample forecast 
load(fullfile(outputFolder, 'qreg_results_inflation.mat'), 'coeffqr_SL'); % qreg 
coeffqreg_SL = permute(coeffqr_SL, [2 1 4 3 5]); % reshape qreg -> quant x var x time x hor x spec 

load(fullfile(outputFolder, 'explanatoryvar_infl.mat'),'explvar_all_spec', 'exovar'); % load explanatory variables acrosss all specifications and dummy covid 
predictors = [ones(size(exovar)), explvar_all_spec, exovar]; % matrix of predictors 
predictors = predictors(168:end,:); % only for the period of interest (2003m12 onward --> first forecast is 2004m1 1 period ahead)

% pre-allocate space for predicted_quantiles_check
predicted_quantiles_check = nan(size(predictors,1),size(coeffqreg_SL,1),size(coeffqreg_SL,4)); % time x qt x hor
contribution_pred = nan(size(predictors,1),size(coeffqreg_SL,1),size(coeffqreg_SL,4), size(predictors,2)); % time x qt x hor x predict

% get the predicted quantiles and the contribution of each predictor
for tt = 1:size(predictors,1)

    for hh = 1:horizons

        pred_q = predictors(tt,:) * squeeze(coeffqreg_SL(:,:,tt,hh))' ; % predicted quantiles
        [sorted_pred, order] = sort(pred_q, 2); % sort along quantile dim

        predicted_quantiles_check(tt,:,hh) = sorted_pred;
        
        for vv = 1:size(predictors,2)
            contri  = predictors(tt,vv) * squeeze(coeffqreg_SL(:,vv,tt,hh))' ; % contribution of each predictor 
            contribution_pred(tt,:,hh,vv) = contri(order);
        end

    end
end

% predicted_quantiles_check == to predicted_quantiles_SL ok
% check that contribution_pred sum to predicted_quantiles_check

%squeeze(contribution_pred(2,:,1,:))
%q5  --> -0.4102    0.5348    0.8834   -0.0141    0.0003   -0.1817         0  =  0.8125
%q10 --> -0.6719    0.5370    1.1149   -0.0956    0.0081   -0.0305         0  =  0.8620
%q15 --> -0.4281    0.5283    1.0688   -0.0679    0.0026   -0.0906         0  =  1.0131
%predicted_quantiles_check(2,:,1); 
%0.8125    0.8620    1.0132  ok
% predicted_quantiles_SL(2,:,1)
% 0.8125    0.8620    1.0132 ok

%%%%%%%%%
% CHART %
%%%%%%%%%

% convert to datetime to construct the date variable for the chart
StartEstDT = datetime(StartEst, 'ConvertFrom','datenum'); % first period at which fcst are estimated
endTDT =  datetime(endT, 'ConvertFrom','datenum'); % last period at which the forecast are performed

% get actual inflation (168 is m1-2004)
actual_var_long_filtered = actual_var_long(idx_estimation:end, :); 
covid_dummy_filtered = covid_dummy_infl(idx_estimation:end, :); 

h_list = [2, 4, 13]; % horizon that we want to plot
quantile_of_interest = [5, 15];  % quantiles we want (1 is 0.05 - 2 is 0.10 and so on..]

for qq = 1:numel(quantile_of_interest)
    
    for ii = 1:numel(h_list)

    q_plot = quantile_of_interest(qq); % quantile of interest
    h_plot = h_list(ii);      % horizon of interest   
    
    % inf_to_plot_evol = actual_var_long_filtered(h_plot:end); % shift inflation forward to match the forecast
    
    % prepare the shaded area 
    dummy_covid_to_plot = covid_dummy_filtered(1:end-h_plot); % shift to match the forecast

    months = (StartEstDT + calmonths(h_plot-1)) : calmonths(1) : endTDT;       % get the quaretrs 

    quantile_q = squeeze(predicted_quantiles_SL(1:size(months,2), q_plot, h_plot, spec_to_use));  % quantile of interest
 
    figure('Position',[100 100 1200 800],'Color','w'); % size and white background
    
    tl = tiledlayout(1,1,'TileSpacing','compact','Padding','compact'); 
    ax = nexttile;

    hold on; 

    % get the quantile of interest 
    contrib_mat = squeeze(contribution_pred(1:size(months,2),q_plot,h_plot,1:7));
    bbar = bar(ax, months, contrib_mat, 'stacked', 'BarWidth',0.6);

    % set nice colormap for 7 predictors
    cmap = lines(7);
    for jj = 1:7
        bbar(jj).FaceColor = cmap(jj,:);
    end   

    plot(ax, months, quantile_q, '-', 'Color',[0 0 0.7], 'LineWidth',1,'DisplayName', sprintf('Quantile q%d', q_plot*5)); % quantile

    yline(0, 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off'); % zer line

    yl = ylim;  % grab them for our patch

    % plot the shaded area
    idx1 = find(dummy_covid_to_plot==1); % find the first occurrence of covid

    if ~isempty(idx1)

        startIdx = idx1(1); % forecast done at 2020q2
        endIdx   = idx1(end); % forecast done at 2022q2
        ms = months(startIdx);
        me = months(endIdx);
        yl = ylim;  % get current y‑limits
        
        % draw exactly one grey patch
        patch([ms me me ms], ...
              [yl(1) yl(1) yl(2) yl(2)], ...
              [0.6 0.6 0.6], 'EdgeColor','none', ...
              'FaceAlpha',0.3, 'HandleVisibility','off');
        uistack(findobj(gca,'Type','patch'),'top'); % change top or bottom
    end

    % plot(months, inf_to_plot_evol, 'k','DisplayName','Outturn','LineWidth',1.25); % actual var
    
    % ylim([0, inf]);
    grid on;

    ticks = months(1):calyears(2):months(end); % every 2 years
    ax.XTick = ticks;
    ax.XAxis.TickLabelFormat = 'yyyy';   % show only the year
    
    % legend above the axes
    lgd = legend(ax, [bbar, findobj(ax,'DisplayName',sprintf('Quantile q%d', q_plot*5))], ...
        [varnames_plot, {sprintf('Quantile q%d', q_plot*5)}], ...
        'Location','northoutside','Orientation','horizontal');
    lgd.Box = 'off';
    set(lgd,'FontSize',18);
    
    set(gca,'FontSize',18);

    hold off;

    end 

end


% end :)
