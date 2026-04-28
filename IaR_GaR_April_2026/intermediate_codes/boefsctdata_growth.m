                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %%%       DATA CLEANING    %%%
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Forecast distribution parameters

% First forecast done in August 2004 --> first forecasted period is q3 2004
% Last forecast done in February 2024 --> last forecasted period is q1 2027
start = datetime('30-Sep-2004', 'InputFormat', 'dd-MMM-yyyy'); % last date with actual data

% Calculate the number of quarters elapsed between start_date and end_fcst
distancetostart = (year(start_date) - year(start)) * 4 + (ceil(month(start_date)/3) - ceil(month(start)/3));
quarters = (year(end_date) - year(start_date)) * 4 + (ceil(month(end_date)/3) - ceil(month(start_date)/3)); 
quarters_fcst = (year(end_fcst) - year(start_date)) * 4 + (ceil(month(end_fcst)/3) - ceil(month(start_date)/3));

% gives the Excel columns. 
col_letter_end = numberToExcelColumn(7+quarters_fcst+distancetostart);
col_letter_start = numberToExcelColumn(7+distancetostart);

% get the rows letter
row_numb_start = (distancetostart)* 10 + 22;
row_numb_end = (distancetostart + quarters +1)* 10 + 22-1;

% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", quarters_fcst +1); 

% Specify sheet and range
opts.DataRange = sprintf('%s%d:%s%d', col_letter_start, row_numb_start, col_letter_end, row_numb_end);
opts.Sheet = 'GDP Forecast';

% Specify column names and types
opts.VariableTypes = repmat("double", 1, quarters_fcst + 1);

% Import the data
growthprojectionparametersmpcS1 = readtable([cd '\Data\kcl_data\gdp_growth_projection_parameters_mpc.xlsx'], opts, "UseExcel", false);

% clear opts 
clear opts

% convert to array 
fcst_growth = table2array(growthprojectionparametersmpcS1);

% Note:
% We focus on forecast based on market interest rate
% we extract mean, mode and sigma (the sqrt of the dispersion)

% col number 
numCols = size(fcst_growth, 2);

% Extract every 10th row for each variable group using direct indexing

% MODE
idx_mode = 1:10:size(fcst_growth, 1);
mrkt_mode_raw = fcst_growth(idx_mode, :)';  
mrkt_mode = reshape(mrkt_mode_raw(~isnan(mrkt_mode_raw)), 13, []);

% MEAN.
idx_mean = 3:10:size(fcst_growth, 1);
mrkt_mean_raw = fcst_growth(idx_mean, :)';
mrkt_mean = reshape(mrkt_mean_raw(~isnan(mrkt_mean_raw)), 13, []);

% UNCERTAINTY (SIGMA)
idx_uncer = 4:10:size(fcst_growth, 1);
mrkt_uncer_raw = fcst_growth(idx_uncer, :)';
mrkt_uncer = reshape(mrkt_uncer_raw(~isnan(mrkt_uncer_raw)), 13, []);

%% final cleaning  

mtestdata = [];

% Loop through each column index
for i = 1:size(mrkt_mode, 2)
    % Append the columns from A, B, and C to the result matrix
    mtestdata = [mtestdata mrkt_mode(:, i) mrkt_mean(:, i) mrkt_uncer(:, i)];
end

clearvars -except mtestdata start_date end_date covid_date ctrynames varnames fullFileName momentlist outputFolder covid_end_date



