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

T = countryData.(ctrynames{1});

% 1) Parse your “yyyymM” strings into datetimes on the 1st of each month
rawDates = T.Date;                          % cell-array of e.g. {'2004m9';'2004m10';…}
yearVec  = str2double(extractBetween(rawDates,1,4));
monthVec = str2double(extractAfter  (rawDates,'m'));
dateVec  = datetime(yearVec, monthVec, 1); % gives e.g. 2024-09-01, 2024-10-01, …

% 2) Filter to the period 2004-09-01 through 2022-09-30
% start_date = datetime(2004,1,1);
% end_date   = datetime(2024,12,1);
keepIdx    = (dateVec >= start_date) & (dateVec <= end_fcst);

T_filtered = T(keepIdx, :);
dateVec    = dateVec(keepIdx);             % now only the months you care about

% 3) (Optional) display these dates as “01-Sep-2004” if you need text
dateVec.Format = 'dd-MMM-yyyy';            % now dateVec(i) displays as '01-Sep-2004'

% pull out the numeric series
arrayvar = T_filtered.(sheetName);   % e.g. 1×N vector of your monthly inflations

% build your lag-matrix exactly as before
n = numel(arrayvar);
interm = NaN(n,n);
for t = 1:n
    interm(1:(n-t+1), t) = arrayvar(t:end);
end

% take the first 37 rows → your 0- to 36-month-ahead actuals
actualvar = interm(1:37, :);

clearvars -except mtestdata start_date end_date covid_date ctrynames varnames fullFileName actualvar momentlist outputFolder
