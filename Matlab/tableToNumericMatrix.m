function data_m = tableToNumericMatrix(tbl)
% tableToNumericMatrix Converts a table to a numeric matrix, handling strings and categoricals.
    nRows = height(tbl);
    nCols = width(tbl);
    data_m = nan(nRows, nCols); % Preallocate with NaN

    for col = 1:nCols
        colData = tbl{:, col};
        if isnumeric(colData)
            data_m(:, col) = colData;
        elseif iscellstr(colData) || isstring(colData)
            % Convert strings to double, empty strings to NaN
            temp = str2double(strtrim(string(colData)));
            data_m(:, col) = temp;
        elseif iscategorical(colData)
            temp = str2double(string(colData));
            data_m(:, col) = temp;
        else
            % For other types, try to convert or leave as NaN
            try
                data_m(:, col) = double(colData);
            catch
                % Leave as NaN
            end
        end
    end
end
