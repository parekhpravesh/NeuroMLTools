function [location, location_U, location_L, cutoff_U, cutoff_L] = detect_outliers(matrix, method, threshold)
% Function to detect outliers within columns of a matrix, given a method
%% Inputs:
% matrix:       matrix to detect outliers in (for each column)
% method:       should be one of:
%                   * 'SD'
%                   * 'IQR'
%                   * 'MAD'
%                   * 'percentile'
% threshold:    number(s) controlling which values are identified as
%               outlier(s) (see Notes)
% 
%% Outputs:
% location:     a logical vector indicating, for each column, the locations 
%               which have been detected as outliers
% location_U:   a logical vector indicating, for each column, the locations
%               which exceed the upper cutoff value
% location_L:   a logical vector indicating, for each column, the locations
%               which are lower than the lower cutoff value
% cutoff_U:     a vector containing upper threshold value for each column
% cutoff_L:     a vector containing lower threshold value for each column
% 
%% Notes:
% Method: SD
% ----------
% Threshold:    number indicating how many standard deviations away
%               from the mean a value has to be to be called as an outlier
% Outlier:      value > (mean + threshold*SD) OR 
%               value < (mean - threshold*SD)
% 
% Method: IQR
% -----------
% Threshold:    three numbers where:
%                 1st number: how many times away from IQR
%                 2nd number: lower percentile 
%                 3rd number: upper percentile 
% Outliers:     values > (upper percentile + threshold*IQR) OR 
%               values < (lower percentile - threshold*IQR) 
% 
% Method: MAD
% -----------
% Threshold:    number indicating how many times away from scaled median
%               absolute deviation (MAD) a number should be to be called an
%               outlier
% Outliers:     values > threshold*scaled MAD OR
%               values < threshold*scaled MAD
% Scaled MAD:   Scaled median absolute deviation is defined as:
%                   c*median(abs(A-median(A)))
%                   c=-1/(sqrt(2)*erfcinv(3/2)) [approximately 1.4826]
% 
% Method: percentile
% ------------------
% Threshold:    two numbers where:
%                   1st number: lower percentile
%                   2nd number: upper percentile
% Outliers:     values > upper threshold (upper percentile) OR 
%               values < lower threshold (lower percentile) OR 
%  
%% References:
% https://www.mathworks.com/help/matlab/ref/isoutlier.html
% https://en.wikipedia.org/wiki/Median_absolute_deviation
% 
%% Defaults:
% method:       'IQR'
% threshold:    3           (for SD)
%               [1.5 25 75] (for IQR)
%               3           (for MAD)
%               [10 90]     (for percentile)
% 
%% Author(s):
% Parekh, Pravesh
% December 23, 2019
% Updated April 22, 2021
% MBIAL

%% Check inputs and assign defaults
% Check matrix
if ~exist('matrix', 'var') || isempty(matrix)
    error('Please provide a vector or matrix to work with');
end

% Check method
if ~exist('method', 'var') || isempty(method)
    method = 'iqr';
else
    method = lower(method);
end

% Check threshold
if ~exist('threshold', 'var') || isempty(threshold)
    if strcmpi(method, 'sd')  || strcmpi(method, 'mad')
        threshold = 3;
    else
        if strcmpi(method, 'iqr')
            threshold = [1.5 25 75];
        else
            if strcmpi(method, 'percentile')
                threshold = [10 90];
            end
        end
    end
else
    % Make sure correct number of threshold values are present
    if (strcmpi(method, 'sd') || strcmpi(method, 'mad'))
        if length(threshold) ~= 1
            error('Only one threshold value required for SD or MAD-based outlier detection');
        end
    else
        if strcmpi(method, 'percentile') 
            if length(threshold) ~= 2
                error('Lower and upper percentile values needed');
            end
        else
            if length(threshold) ~= 3
                error('Times away from IQR, and lower and upper percentile values needed');
            end
        end
    end
end

% Check if Statistics and Machine Learning Toolbox is present
if license('test','Statistics_toolbox')
    useStats = true;
else
    useStats = false;
end

%% Set upper and lower cutoff
switch(method)
    case 'sd'
        all_mean    = mean(matrix, 1);
        all_SD      = std(matrix, [], 1);
        cutoff_U    = all_mean + threshold * all_SD;
        cutoff_L    = all_mean - threshold * all_SD;

    case 'iqr'
        if useStats
            prctile_L   = prctile(matrix, threshold(2));
            prctile_U   = prctile(matrix, threshold(3));
            all_IQR     = iqr(matrix);
            cutoff_U    = prctile_U + threshold(1) * all_IQR;
            cutoff_L    = prctile_L - threshold(1) * all_IQR;
        else
            [~, ~, ~, prctLU, IQRs] = calc_percentiles(matrix, [threshold(2), threshold(3)]);
            cutoff_L                = (prctLU(:,1) - threshold(1) * IQRs)';
            cutoff_U                = (prctLU(:,2) + threshold(1) * IQRs)';
        end
        
    case 'mad'
        c = -1/(sqrt(2)*erfcinv(3/2));
        if useStats
            all_median  = median(matrix, 1);
            tmp         = bsxfun(@minus, matrix, all_median);
            MAD         = median(abs(tmp));
        else
            [~, all_median] = calc_percentiles(matrix);
            all_median      = all_median';
            tmp             = bsxfun(@minus, matrix, all_median);
            [~, MAD]        = calc_percentiles(abs(tmp));
            MAD             = MAD';
        end
        cutoff_U    = all_median + (threshold * c * MAD);
        cutoff_L    = all_median - (threshold * c * MAD);
        
    case 'percentile'
        if useStats
            cutoff_L    = prctile(matrix, threshold(1));
            cutoff_U    = prctile(matrix, threshold(2));
        else
            [~, ~, ~, cutoff] = calc_percentiles(matrix, [threshold(1), threshold(2)]);
            cutoff_L          = cutoff(:,1)';
            cutoff_U          = cutoff(:,2)';
        end
        
    otherwise
        error(['Unknown outlier detection method specified: ', method]);
end

%% Mark outliers
location_U  = bsxfun(@gt, matrix, cutoff_U);
location_L  = bsxfun(@lt, matrix, cutoff_L);
location    = location_U | location_L;