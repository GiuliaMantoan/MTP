                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    %%%       DATA CLEANING    %%%
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% empty strucutre to store country data  
countryData = struct();
countryData.(ctrynames{1}) = table();

sheetName = varnames{1}; % select the variable
try
    T = readtable(fullFileName, 'Sheet', sheetName);
catch ME
    error('Error reading the sheet "%s": %s', sheetName, ME.message);
end
    
% extract data

country = ctrynames{1};

% Find the index of the country column
colIdx = find(strcmpi(T.Properties.VariableNames, ctrynames{1}));

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

%% clean the data

%% clean the data (now with monthly, quarter-end sub-sampling)

T = countryData.(ctrynames{1});

% parse your “2004m9”-style strings into true datetimes
nObs = height(T);
dateVec = NaT(nObs,1);
for k = 1:nObs
    dStr = char(T.Date{k});          % e.g. '2004m9' or '2004m10'
    yearNum    = str2double(dStr(1:4));
    monthNum   = str2double(dStr(6:end));
    % use day=1 since we only care about month-ends by month number
    dateVec(k) = datetime(yearNum, monthNum, 1);
end

% define your window
start_date = datetime(2004, 9, 1);
end_date   = datetime(2022, 9, 30);

% keep only quarter-ends (Mar, Jun, Sep, Dec) in that window
isQE = ismember(month(dateVec), [3 6 9 12]);
inWindow = (dateVec >= start_date) & (dateVec <= end_date);
keepIdx = isQE & inWindow;

T_filtered = T(keepIdx, :);
filteredDates = dateVec(keepIdx);

% now pull out the numeric data
arrayvar = T_filtered.(sheetName);   % one-column of values

% build your lag-matrix exactly as before
n = numel(arrayvar);
interm = NaN(n,n);
for i = 1:n
    interm(1:(n-i+1), i) = arrayvar(i:end);
end

% take the first 13 lags → your "actualvar"
actualvar_wlrt = interm(1:13, :);

clearvars -except mtestdata start_date end_date covid_date ctrynames varnames fullFileName actualvar_wlrt actualvar momentlist outputFolder modelfcst covid_qtr covid_idx
