% Authors: David Aikman, Rhys Bidder, Simone Maso, Aditya Mori 

%% Data

% covid treatment 
% reminder for Adi: if you are using this .m file and gives you an error
% here. it is due to this covid dummy. 
% if you want to be coherent with what you had before just set covid = 0. 

if covid == 1
    varnames = [ varnames , dummyvarname ]; % placeholder for dummyvarname
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

% %%%%% UNCOMMENT FOR OUTERJOIN %%%%%%%%%
% 
%         else
%     % Build a tiny table for the new series
%     newT = table(dates, data, 'VariableNames', {'Date', sheetName});
% 
%     % Outer-join on Date so rows align and missing months become NaN
%     countryData.(country) = outerjoin(countryData.(country), newT, ...
%         'Keys', 'Date', 'MergeKeys', true, 'Type', 'full');
% 
%     % (No suffix clean-up needed because we used MergeKeys=true and
%     % the new variable keeps its name "sheetName")
%             % --- Ensure rows are in true time order and unique by month ---
% tmp = countryData.(country);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ---------------------------------------------------------------
% % Create a numeric month key from 'YYYYmM' / 'YYYYmMM',
% % then sort the table by that key and remove duplicate months.
% % This ensures rows are in true chronological order with one row per month.
% % ---------------------------------------------------------------
% 
% % 1) Build a numeric date key (e.g., 1990m2 → datenum(1990,2,1))
% dateNumeric_tmp = zeros(height(tmp),1);        % preallocate (one key per row)
% for kk = 1:height(tmp)
%     dStr = char(tmp.Date{kk});                 % raw date string, e.g. '1990m2' or '1990m12'
%     yr   = str2double(dStr(1:4));              % '1990' → 1990
%     mn   = str2double(dStr(6:end));            % '2' or '12' → 2 or 12
%     dateNumeric_tmp(kk) = datenum(yr, mn, 1);  % first day of that month (serial number)
% end
% tmp.dateNumeric = dateNumeric_tmp;             % attach helper key to the table
% 
% % 2) Sort rows by time, then drop duplicate months (keep the first)
% tmp = sortrows(tmp, 'dateNumeric');            % ascending chronological order
% [~, ia] = unique(tmp.dateNumeric, 'stable');   % indices of first occurrence per month
% tmp = tmp(ia, :);                              % keep only the first row for each month
% 
% % 3) Clean up the helper column now that sorting/dedup is done
% tmp.dateNumeric = [];
% 
% % 4) Write the cleaned table back into the countryData struct
% countryData.(country) = tmp;
% 
% % 5) Tidy temporary variables
% clear tmp dateNumeric_tmp
% end
%       end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
        yearNum  = str2double(dStr(1:4));    % "1970" → 1970
        monthNum = str2double(dStr(6:end));  % skip the 'm', take "6" or "12"
        dateNumeric(k) = datenum(yearNum, monthNum, 1);
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
    yearNum  = str2double(dStr(1:4));    % "1970" → 1970
    monthNum = str2double(dStr(6:end));  % skip the 'm', take "6" or "12"
    dateNumeric(k) = datenum(yearNum, monthNum, 1);

end


% clean the environment  
% clearvars -except countryData ctrynames ctryweightsname endT lags startT weightedvarnames varnames fullFileName weightsFile outputFolder dateNumeric StartEst onlyuk quantilelevels model_selection modellist horizons spec combo_specifications
clearvars  c colIdx d data dateNumeri dates idx j k laggedData lagglobalvar m quarterNum newVarName yearNum v T sheetName lag i dStr