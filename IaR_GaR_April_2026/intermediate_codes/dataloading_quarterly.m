% Authors: David Aikman, Rhys Bidder, Simone Maso, Aditya Mori 

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
clearvars  c colIdx d data dateNumeri dates idx j k laggedData lagglobalvar m quarterNum newVarName yearNum v T sheetName lag i dStr