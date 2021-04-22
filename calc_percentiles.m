function [percentile25, percentile50, percentile75, percentileot, IQRs] = calc_percentiles(matrix, percentiles)
% Function that calculates 25th, 50th, and 75th percentiles, along with
% other desired percentiles for each column of a given matrix
%% Inputs:
% matrix:           [n x p] matrix for which to calculate percentiles for 
%                   each column p
% percentiles:      [k x 1] vector indicating which k percentiles are to be 
%                   calculated for each column p
% 
%% Outputs:
% percentile25:     [p x 1] 25th percentile for each column
% percentile50:     [p x 1] 50th percentile or median for each column 
% percentile75:     [p x 1] 75th percentile for each column
% percentileot:     [p x k] matrix of entries corresponding to the k 
%                   percentiles calculated on each column p
% IQRs:             [p x 1] difference between 75th and 25th percentiles
%                   for each column
% 
%% Notes:
% The percentile detection method is based on the description provided in
% the Statistics and Machine Learning Toolbox:
% 
% For each column p in matrix:
%   * Find the number of non-NaN entries n
%   * Sort the elements in ascending order
%   * Assign these percentiles to the sorted elements:
%       (100 * (0.5:1:n))/n
%   * Use linear interpolation for percentile values between these ranges
%   * For percentile values beyond this range, assign minimum or maximum of
%     the values in p
% 
% Linear interpolation is performed as described in the documentation for
% prctile the Statistics and Machine Learning Toolbox:
% 
% For each column p in matrix:
%   For each percentile value:
%       * Let y1 be the first percentile value which is lower than the 
%         required percentile; let yy1 be the percentile at this point
%       * Let y2 be the first percentile value which is higher than the 
%         required percentile; let yy2 be the percentile at this point
%       * Let yy be the required percentile value to be calculated
%       * New value = y1 + ((yy - yy1)/(yy2 - yy1)) * (y2 - y1)
% 
% The 50th percentile is calculated as the mean of the middle value:
% For each column p in matrix:
%   * If n, the number of non-NaN entries are even
%       * median = average of n/2 and (n/2 + 1) values
%   * If n, the number of non-NaN entries are odd
%       * median = value at (n+1)/2
% 
%% References:
% MATLAB documentation for calculating percentiles:
% https://www.mathworks.com/help/stats/prctile.html
% 
%% Defaults:
% percentiles:  []
% 
%% Author(s):
% Parekh, Pravesh
% April 21, 2021
% MBIAL

%% Check inputs
% Check matrix
if ~exist('matrix', 'var') || isempty(matrix)
    error('Please provide a matrix to calculate percentiles on');
else
    num_vars = size(matrix,2);
end

% Check percentiles
if ~exist('percentiles', 'var') || isempty(percentiles)
    percentiles  = [];
end

%% Initialize
percentile25 = zeros(num_vars, 1);
percentile50 = zeros(num_vars, 1);
percentile75 = zeros(num_vars, 1);
percentileot = zeros(num_vars, length(percentiles));
IQRs         = zeros(num_vars, 1);

%% Calculate percentiles
for vars = 1:num_vars
    % Locate NaNs
    loc_nan  = isnan(matrix(:,vars));
    
    % Number of non-NaN values
    num_vals = sum(~loc_nan);
    
    % Arrange values in ascending order
    tmp_vals = sort(matrix(~loc_nan,vars), 'ascend');
    
    % Assign a percentile to these values
    ass_prct = (100 * (0.5:1:num_vals))/num_vals;
    
    % Find if 25th or 75th percentiles already exist
    loc25 = ass_prct == 25;
    loc75 = ass_prct == 75;
    if sum(loc25) == 1
        percentile25(vars,1) = tmp_vals(loc25);
    else
        % Find first upper value just above 25th percentile; lower will be
        % value just previous to this as the data is sorted
        loc_upper = find(ass_prct > 25, 1);
        loc_lower = loc_upper - 1;
        
        % Interpolate
        percentile25(vars,1) = do_interp(tmp_vals(loc_lower), tmp_vals(loc_upper), ...
                                         ass_prct(loc_lower), ass_prct(loc_upper), 25);
                                         
    end
    
    if sum(loc75) == 1
       percentile75(vars,1) = tmp_vals(loc75);
    else
        % Find first upper value just above 75th percentile; lower will be
        % value just previous to this as the data is sorted
        loc_upper = find(ass_prct > 75, 1);
        loc_lower = loc_upper - 1;
        
        % Interpolate
        percentile75(vars,1) = do_interp(tmp_vals(loc_lower), tmp_vals(loc_upper), ...
                                         ass_prct(loc_lower), ass_prct(loc_upper), 75);
    end 
    
    % Calculate 50th percentile
    if rem(num_vals,2)
        % Is an odd number
        percentile50(vars,1) = tmp_vals((num_vals+1)/2);
    else
        % Is an even number
        percentile50(vars,1) = (tmp_vals(num_vals/2) + tmp_vals(num_vals/2 + 1))/2;
    end
    
    % Calculate IQR
    IQRs(vars,1) = percentile75(vars,1) - percentile25(vars,1);
    
    % Calculate remaining percentiles
    for prc = 1:length(percentiles)
        
        % Check if the value is already calculated
        loc = ass_prct == percentiles(prc);
        if sum(loc) == 1
            percentileot(vars, prc) = tmp_vals(loc);
        else
            loc_upper = find(ass_prct > percentiles(prc),1);
            loc_lower = loc_upper - 1;
            
            % If upper value does not exist, then assign the maximum value
            % of this particular column as this percentile
            if isempty(loc_upper)
                percentileot(vars,prc) = max(tmp_vals);
            else
                % If lower value does not exist, then assign the minimum
                % value of this particular column as this percentile
                if isempty(loc_lower) || loc_lower == 0
                    percentileot(vars,prc) = min(tmp_vals);
                else
                    % Interpolate
                    percentileot(vars,prc) = do_interp(tmp_vals(loc_lower), tmp_vals(loc_upper), ...
                                                       ass_prct(loc_lower), ass_prct(loc_upper), percentiles(prc));
                end
            end
        end
    end
end

function value = do_interp(y1, y2, yy1, yy2, yy)
value = y1 + ((yy - yy1)/(yy2 - yy1)) * (y2 - y1);