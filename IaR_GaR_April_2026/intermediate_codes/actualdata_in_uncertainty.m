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

% Convert the 'Date' column to datetime using quarterly parsing
quarterStrs = T.Date;

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


% Logical index for dates between start_date and end_date
keepIdx = dateNumeric >= datenum(start_date) & dateNumeric <= datenum(end_date);

% Filter table
T_filtered = T(keepIdx', :);

% get array variable
arrayvar = T_filtered.(sheetName);

% Number of observations
n = length(arrayvar);

% Initialize a matrix to hold the columns
interm = NaN(n, n);

% Fill the matrix with the required observations
for i = 1:n
    interm(1:(n-i+1), i) = arrayvar(i:end);
end

actualvar = interm(1:13,:);

%clearvars -except mtestdata start_date end_date covid_date ctrynames varnames fullFileName actualvar momentlist outputFolder
