% Authors: David Aikman, Rhys Bidder, Simone Maso, Aditya Mori 

% close all; clear; clc;

%% Settings
 
%set(0,'defaultAxesFontName', 'Times'); % font for chart axis 
%set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1) % line style 

%rng('default'); % defaul random number generateor 
%rng(0);     

%cd 'C:\Users\k2370179\Dropbox\BoE-KCL Macro Forecasting\' % chage here the cd 
%addpath('Data\') % data 
% addpath('functions') % functions % addpath('intermediate_code') % functions 
%addpath('Outputs\') % functions 

% file names
%fullFileName = fullfile('Data', 'GaRDataRaw.xlsx'); 
%weightsFile = fullfile('Data', 'ExportWeights.xlsx');

%varnames = {'g4Infl', 'InfExp', 'ygap'}; % variables you want to import 
%lags = [1 0 0 2 0]; % this create lagged variables and add them to the dataset (each element is the number of lag) the last two refers to weighted_exp and gr4_weighted_exp  
%ctrynames = {'UK', 'USA'}; % country you want to import 
%weightedvarname = 'PX'; % weighted var
%ctryweightsname = {'AUS','CAN','DNK','JAP','NOR','SWE','SWI','UK','USA','KOREA','EA'}; % country to create world export prices 
%startT = datenum(1997,03,31); % sample start 
%endT   = datenum(2022,09,30); % sample end 

%% Data

% covid treatment 
% reminder for Adi: if you are using this .m file and gives you an error
% here. it is due to this covid dummy. 
% if you want to be coherent with what you had before just set covid = 0. 

if covid == 1
    varnames = [ varnames , {'covid'} ];
    lags = [lagctryvar 0 lagglobalvar];

end

% empty strucutre to store country data  
countryData = struct();
for i = 1:length(ctrynames)
    countryData.(ctrynames{i}) = table();
end

for v = 1:length(varnames)

    sheetName = varnames{v}; % select the variable

    try
        T = readtable(fullFileName, 'Sheet', sheetName);
    catch ME
        error('Error reading the sheet "%s": %s', sheetName, ME.message);
    end
    
    % Loop over each country to extract data
    for i = 1:length(ctrynames)

        country = ctrynames{i};

        % Find the index of the country column
        colIdx = find(strcmpi(T.Properties.VariableNames, country));

        if isempty(colIdx)
            error('Country "%s" not found in the sheet "%s".', country, sheetName);
        end
        
        % Extract the date and the country data
        dates = T{:,1};    
        data  = T{:,colIdx};
        
        % If this is the first variable being processed, create a new table
        if isempty(countryData.(country))
            countryData.(country) = table(dates, data, 'VariableNames', {'Date', sheetName});
        else
            % Otherwise, add a new variable column to the existing table
            countryData.(country).(sheetName) = data;
        end
    end
end

%% Generate weighted export

%{

% Empty structure to store weighting data  
weightData = struct();

for i = 1:length(ctrynames)
    weightData.(ctrynames{i}) = table();
end

% Loop over each country in ctry
for i = 1:length(ctrynames)

    % Construct the sheet name by appending "_X" to the weighting country code
    sheetName = [ctrynames{i}, '_X'];
    
    % Read the weighting data from the specified sheet
    try
        W = readtable(weightsFile, 'Sheet', sheetName);
    catch ME
        error('Error reading the sheet "%s": %s', sheetName, ME.message);
    end
    
    % Restrict the data to only the valid countries from the country_weights list.
    selectedCols = ismember(W.Properties.VariableNames, ctryweightsname);
    
    if ~any(selectedCols)
        error('None of the columns in the sheet "%s" match the expected countries: %s.', ...
            sheetName, strjoin(ctryweightsname, ', '));
    end

    % Extract only the selected columns and store them in the structure
    weightData.(ctrynames{i}) = [W(:,2) W(:, selectedCols)];

end

for vv = 1:length(weightedvarnames)

    weightedvarname = weightedvarnames{vv};

    % get the weighted data - in our case export prices
    try
        WW = readtable(fullFileName, 'Sheet', weightedvarname);
    catch ME
        error('Error reading the sheet "%s" from GaRDataRaw.xlsx: %s', weightedvarname, ME.message);
    end

    % Restrict the WW table to only the columns corresponding to the valid weighting countries
    selectedColsWW = ismember(WW.Properties.VariableNames, ctryweightsname);
    if ~any(selectedColsWW)
        error('None of the columns in the sheet "%s" of GaRDataRaw.xlsx match the expected countries: %s.', ...
            weightedvarname, strjoin(ctryweightsname, ', '));
    end

    % Create a new table with the date and selected weighting columns
    WW_filtered = [WW(:,1), WW(:, selectedColsWW)];

    % construct the weighted var for each country and add it to the country
    % data tables

    for i = 1:length(ctrynames)
    
        currentCountry = ctrynames{i};
        
        % Get the weighting table for the current country (weightData already contains the date column)
        W_table = weightData.(currentCountry);
        
        % Extract the date vector from weightData (assumed to be the first column)
        weights_dates = W_table{:,1};
        
        % Get the names of the weight columns 
        weightCols = W_table.Properties.VariableNames(2:end);
    
        % Reorder weights columns to match the order in country_weights
        [~, idx_order] = ismember(ctryweightsname, weightCols);
        if any(idx_order == 0)
            error('Some weighting columns are not found in weightData for %s.', currentCountry);
        end
        weights_matrix = W_table{:, 1+idx_order};
    
        % replace NaN with zeros
        weights_matrix(isnan(weights_matrix)) = 0;
    
        % Extract the weighted variable from WW_filtered (first column is date)
        weightedvar_dates = WW_filtered{:,1};
        weightedvarCols = WW_filtered.Properties.VariableNames(2:end);
        % Reorder export price columns to match the order in country_weights
        [~, idx_order2] = ismember(ctryweightsname, weightedvarCols);
        if any(idx_order2 == 0)
            error('Some export price columns are not found in WW_filtered.');
        end
        weightedvar_matrix = WW_filtered{:, 1+idx_order2};
        
        % Check that the dates match between weighting data and export prices
        if ~isequal(weights_dates, weightedvar_dates)
            error('Date mismatch between weightData and WW_filtered for country %s.', currentCountry);
        end
        
        % Compute the weighted export price for each row:
        weightedvar = sum(weights_matrix .* weightedvar_matrix, 2) ./ sum(weights_matrix, 2);
        
        % Add the computed weighted export as a new column to the corresponding countryData table.
        countryData.(currentCountry).(weightedvarname) = weightedvar;

        if ygrowthrate == 1
    
        % Compute year-on-year growth rate for weighted_exp.
        grweighted_var_4 = NaN(size(weightedvar));  % Initialize with NaNs
        grweighted_var_4(5:end) = 100 * (weightedvar(5:end) - weightedvar(1:end-4)) ./ weightedvar(1:end-4);
        
        % Add the growth rate column to the corresponding countryData table.
        growthVarName = ['g4' weightedvarname];
        countryData.(currentCountry).(growthVarName) = grweighted_var_4;

        end
    
    end
  
end

% double check UK 1995q1: 83.5485 OK

%}

%% Variables selection and transformation and time selection 

% Loop over each country's table
for c = 1:length(ctrynames)

    % Get the current table
    T = countryData.(ctrynames{c});
    
    % Loop over each variable column after the date
    for j = 2:width(T)

        % Determine how many lags to create for this variable using the lag vector
        for lag = 1:lags(j-1)

            % Create the lagged variable; pad the beginning with NaNs
            laggedData = [nan(lag, 1); T{1:end-lag, j}];

            % Generate a new variable name
            newVarName = sprintf('l%d%s', lag, T.Properties.VariableNames{j});

            % Add the new column to the table

            T.(newVarName) = laggedData;
        end
    end
    
    % --- Now filter the data by date ---
    dates = T{:,1};
    dateNumeric = zeros(length(dates),1);

    for k = 1:length(dates)
        % Convert the date string to a character array
        dStr = char(dates{k});
        yearNum = str2double(dStr(1:4));
        quarterNum = str2double(dStr(end));
        % Map the quarter to a month and day (quarter-end dates)
        switch quarterNum
            case 1, m = 3;  d = 31;
            case 2, m = 6;  d = 30;
            case 3, m = 9;  d = 30;
            case 4, m = 12; d = 31;
            otherwise, error('Unexpected quarter in date string %s', dStr);
        end
        dateNumeric(k) = datenum(yearNum, m, d);
    end
    
    % Find rows that fall within the desired date range
    idx = (dateNumeric >= startT) & (dateNumeric <= endT);
    
    % Subset the table to keep only the selected rows
    T = T(idx, :);
    
    % Save the updated table back to the struct
    countryData.(ctrynames{c}) = T;
end

% get the index for estimation start
clear dateNumeric
dates = countryData.(ctrynames{1}).Date;

for k = 1:length(dates)
    % Convert the date string to a character array
    dStr = char(dates{k});
    yearNum = str2double(dStr(1:4));
    quarterNum = str2double(dStr(end));
    % Map the quarter to a month and day (quarter-end dates)
    switch quarterNum
            case 1, m = 3;  d = 31;
            case 2, m = 6;  d = 30;
            case 3, m = 9;  d = 30;
            case 4, m = 12; d = 31;
            otherwise, error('Unexpected quarter in date string %s', dStr);
    end
        dateNumeric(k) = datenum(yearNum, m, d);

end


% clean the environment  
% clearvars -except countryData ctrynames ctryweightsname endT lags startT weightedvarnames varnames fullFileName weightsFile outputFolder dateNumeric StartEst onlyuk quantilelevels model_selection modellist horizons spec combo_specifications
clearvars  c colIdx d data dateNumeri dates idx j k laggedData lagglobalvar m quarterNum newVarName yearNum v T sheetName lag i dStr