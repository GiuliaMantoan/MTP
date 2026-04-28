function colLetter = numberToExcelColumn(n)
% numberToExcelColumn converts a positive integer to its Excel column letter.
%
%   colLetter = numberToExcelColumn(n) returns a character vector
%   representing the Excel column letters corresponding to the
%   positive integer n.
%
%   Example:
%       numberToExcelColumn(100) returns 'CV'
%
%   Input:
%       n - Positive integer representing the column number.
%
%   Output:
%       colLetter - A character vector representing the Excel column letters.

    % Validate input
    if n < 1 || floor(n) ~= n
        error('Input must be a positive integer.');
    end

    % Initialize the output as an empty character array.
    colLetter = '';
    
    % Convert the number to letters.
    while n > 0
        remainder = mod(n - 1, 26);
        colLetter = [char(65 + remainder) colLetter];  %#ok<AGROW>
        n = floor((n - remainder - 1) / 26);
    end
end