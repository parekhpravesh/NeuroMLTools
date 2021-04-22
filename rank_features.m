function [rankedList, ranksOrdered] = rank_features(data, classes, rank_method, ...
                                                    std_method, out_method,     ...
                                                    out_thresh, out_handle)
% Function to calculate rank of features given a method
%% Inputs:
% data:         a matrix with rows having samples and columns as features
% classes:      a vector with 0's and 1's corresponding to sample classes
% rank_method:  feature ranking method; should be one of:
%                   * 'tstats'
%                   * 'wilcoxon'
%                   * 'bhattacharyya'
%                   * 'relieff'
%                   * 'mrmr'
%                   * 'dmean'
%                   * 'dmedian'
%                   * 'dstd'
%                   * 'correlation'
%                   * 'cohensd'
%                   * 'hedgesg'
%                   * 'glassd'
%                   * 'silhouetteman'
%                   * 'silhouetteeuc'
% std_method:   standardization method; should be one of:
%                   * 'rescale'
%                   * 'mean'
%                   * 'std'
%                   * 'none'
% out_method:   method for dealing with outliers; should be one of:
%                   * 'SD'
%                   * 'IQR'
%                   * 'MAD'
%                   * 'percentile'
%                   * 'none'
% out_thresh:   if out_method is specified (except 'none'), one (SD or
%               MAD), two (percentile), or three (IQR) values for detecting
%               outliers; see detect_outliers for documentation
% out_handle:   handling of outliers; should be one of:
%                   * 'winsorize'
%                   * 'trim'
% 
%% Outputs:
% rankedList:   a vector arranged in ascending order of ranks i.e. 
%               feat_ranks(1) will be the index of the feature which was 
%               ranked 1
% ranksOrdered: a vector corresponding to rank of each feature (in the
%               order of column numbers of data i.e. feat_ordered(1) will 
%               be the rank of feature 1)
% 
%% Notes:
% Ranks are calculated as follows:
%   1. Find outliers, if required (requires detect_outliers)
%   2. Handle outliers (only if outlier detection is true)
%   3. Standardize feature, if required (requires feature_scaling)
%   4. Rank and sort based on rank_method (see below)
% 
% tstats:           rank by two sided absolute value of T Statistics
%                   assuming unequal variance (i.e. Welch test)
% 
% wilcoxon:         rank by absolute value of U statistics, calculated as:
%                   U = ranksum - (nx*(nx+1))/2, where nx=sample size of 
%                   class = 0; use approximate method for ranksum
% 
% bhattacharyya:	rank by absolute value of univariate Bhattacharyya 
%                   distance; complex values are set to NaN
% 
% relieff:          rank using ReliefF with k=10
% 
% mrmr:             rank using mRMR method (requires mRMR_Spearman or 
%                   mRMR_Spearman_alt function)
% 
% dmean:            rank based on absolute difference in the means
% 
% dmedian:          rank based on absolute difference in the median values 
% 
% dstd:             rank based on absolute difference in the standard
%                   deviations of features
% 
% correlation:      rank based on the absolute value of correlation between 
%                   each feature and the class variable
% 
% cohensd:          rank based on the absolute value of Cohen's d
% 
% hedgesg:          rank based on the absolute value of Hedges' g
% 
% glassd:           rank based on the absolute value of Glass's delta
% 
% silhouetteman:    rank based on the overall mean value of Silhouette
%                   values using Manhattan distance as distance metric
% 
% silhouetteeuc:    rank based on the overall mean value of silhouette
%                   values using Euclidean distance as distance metric
% 
% Note that calculating correlation between X and y variables with the
% assumption that y is a member of {0,1} is the same as calculating the
% point-Biserial correlation:
% term1              = (mean(data(classes==1,:)) - mean(data(classes==0,:)))./std(data);
% term2              = sqrt((sum(classes==1)*sum(classes==0))/(length(classes)*(length(classes)-1)));
% pointbiserial_corr = term1 * term2;
% 
% For details of outlier detection, see detect_outliers
% 
% Outlier handling:
% 'winsorize':      settings used for detecting the outliers are also used
%                   for replacing the outlier values
% 'trim':           outliers are excluded from the ranking procedure
% 
% Further notes on winsorization:
% SD:               values > and < (mean +/- threshold*SD) are set equal to
%                   (mean +/- threshold*SD)
% IQR:              values > and < (upper/lower percentile +/- threshold*IQR)
%                   are set equal to (upper/lower percentile +/- threshold*IQR)
% MAD:              values > and < median +/- threshold*scaled MAD are set
%                   equal to median +/- threshold*scaled MAD
% percentile:       values > and < upper/lower percentile are set equal to 
%                   upper/lower percentile
% 
% Further notes on trim:
% Outlier values are set to NaN and the rest of the calculation is
% performed by using the 'omitnan' flag. However, mRMR_Spearman needs to be
% edited to account for NaN values
% 
% Modifications required for mRMR_Spearman:
% 1) Change line 61 to:
%       CI = abs(corr(X, y, 'type', 'Spearman', 'rows', 'complete'));
% 
% 2) Before the main loop (line 70 onwards), add:
%     % Calculate self correlation
%     all_redundancy = corr(X, 'type', 'Spearman', 'rows', 'complete');
% 
% 3) Change line 79 to:
%     CI_array(idx_pointer(i), last_feature) = abs(all_redundancy(features(last_feature), idx_pointer(i))); 
% 
%% References:
% Calculation of U statistics:
% https://in.mathworks.com/help/stats/ranksum.html
% 
% Calculation of mRMR:
% A. Tsanas, M.A. Little, P.E. McSharry: "A methodology for the analysis of
% medical data", Chapter 7 in Handbook of Systems and Complexity in Health, 
% pp. 113-125, Eds. J.P. Sturmberg, and C.M. Martin, Springer, 2013
% https://www.darth-group.com/software
% 
% Calculation of effect sizes:
% Daniël Lakens. 2013. Calculating and reporting effect sizes to facilitate 
% cumulative science: a practical primer for t-tests and ANOVAs. Frontiers
% in Psychology. https://doi.org/10.3389/fpsyg.2013.00863
% 
%% Defaults:
% rank_method:  'relieff'
% std_method:   'std'
% out_method:   'MAD'
% out_handle:   'trim'
% 
%% Author(s)
% Parekh, Pravesh
% January 11, 2020
% MBIAL

%% Parse inputs and assign defaults
% Check data
if ~exist('data', 'var') || isempty(data)
    error('Please provide a data matrix to work with');
else
    num_features = size(data,2);
    num_samples  = size(data,1);
end

% Check classes
if ~exist('classes', 'var') || isempty(classes)
    error('Please provide a vector of class labels having 0s and 1s');
else
    if length(classes) ~= num_samples
        error('Mismatch between number of samples and number of class labels');
    else
        uq_vals = unique(classes);
        if length(uq_vals) ~= 2
            error('Detected more than two classes; only two classes supported');
        else
            if sum(ismember(uq_vals, [0 1])) ~= 2
                error('Classes should be a vector of 0s and 1s');
            end
        end
    end
end

% Check ranking method
if ~exist('rank_method', 'var') || isempty(rank_method)
    rank_method = 'relieff';
else
    rank_method = lower(rank_method);
end

% Check whether to standardize or not
if ~exist('std_method', 'var') || isempty(std_method)
    to_std     = true;
    std_method = 'std';
else
    if strcmpi(std_method, 'none')
        to_std = false;
    else
        to_std = true;
        std_method = lower(std_method);
    end
end

% Check whether outliers should be removed
if ~exist('out_method', 'var') || isempty(out_method)
    to_out      = true;
    out_method  = 'MAD';
else
    if strcmpi(out_method, 'none')
        to_out = false;
    else
        to_out = true;
    end
end

% Check out_thresh
if ~exist('out_thresh', 'var')
    out_thresh = '';
end

% Check out_handle
if ~exist('out_handle', 'var') || isempty(out_handle)
    out_handle = 'trim';
else
    out_handle = lower(out_handle);
    if ~ismember(out_handle, {'winsorize'; 'trim'})
        error(['Incorrect out_handle specified: ', out_handle]);
    end
end

%% Detect and handle outliers, if required
if to_out
    
    % Detect outliers
    [outliers, outlier_U, outlier_L, cutoff_U, cutoff_L] = detect_outliers(data, out_method, out_thresh);
    
    % Handle outliers
    if strcmpi(out_handle, 'winsorize')
        % Winsorize each feature
        for feat = 1:num_features
            data(outlier_U(:,feat), feat) = cutoff_U(feat);
            data(outlier_L(:,feat), feat) = cutoff_L(feat);
        end
    else
        % Trim features by setting to NaN
        data(outliers) = NaN;
    end
end

%% Standardize data, if required
if to_std
    data = feature_scaling(data, std_method);
end

%% Rank features
switch rank_method
    case 'tstats'
        [~, ~, ~, stats]  = ttest2(data(classes==0,:), data(classes==1,:), 'Vartype', 'unequal', 'Tail', 'both');
        to_rank           = stats.tstat;
        [~, rankedList]   = sort(abs(to_rank'), 'descend');
        [~, ranksOrdered] = sort(rankedList,    'ascend');
        
    case 'wilcoxon'
        numx     = sum(classes==0);
        to_rank  = zeros(num_features,1);
        for feat = 1:num_features
            [~, ~, stats] = ranksum(data(classes==0,feat), data(classes==1,feat),  ...
                                   'method', 'approximate', 'tail', 'both');
            to_rank(feat,1) = stats.ranksum - (numx*(numx+1))/2;
        end
        [~, rankedList]   = sort(abs(to_rank), 'descend');
        [~, ranksOrdered] = sort(rankedList,   'ascend');
        
    case 'bhattacharyya'
        % Define matrices
        mat1 = data(classes==0,:);
        mat2 = data(classes==1,:);
        
        % Get all mean values
        mean1       = mean(mat1, 'omitnan');
        mean2       = mean(mat2, 'omitnan');
        diff_mean   = mean1 - mean2;
        
        % Get all covariance values
        cov1     = cov(mat1, 'partialrows');
        cov2     = cov(mat2, 'partialrows');
        avg_cov  = (cov1 + cov2)/2;
        to_rank = zeros(num_features,1);
        
        for feat = 1:num_features
            try
                C1      = cov1(feat, feat);
                C2      = cov2(feat, feat);
                C       = avg_cov(feat, feat);
                term1   = 1/8 * diff_mean(feat) * (C \ diff_mean(feat)');
                term2   = 1/2 * log(det(C)/sqrt(det(C1)*det(C2)));
                tmp     = term1 + term2;
                if ~isreal(tmp)
                    to_rank(feat,1) = NaN;
                else
                    to_rank(feat,1) = term1 + term2;
                end
            catch
                to_rank(feat,1) = NaN;
            end
        end
        [~, rankedList]   = sort(abs(to_rank), 'descend');
        [~, ranksOrdered] = sort(rankedList,   'ascend');
        
    case 'relieff'
        [rankedList, to_rank] = relieff(data, classes, 10, 'method', 'classification');
        [~, b]                = sort(to_rank, 'descend');
        [~, ranksOrdered]     = sort(b,       'ascend');
        rankedList            = rankedList';
        ranksOrdered          = ranksOrdered';
        
    case 'mrmr'
        try
            rankedList    = mRMR_Spearman_alt(data, classes, num_features);
        catch
            warning('Unable to find mRMR_Spearman_alt; trying mRMR_Spearman');
            rankedList    = mRMR_Spearman(data, classes, num_features);
        end
        rankedList        = rankedList';
        [~, ranksOrdered] = sort(rankedList, 'ascend');
        
    case 'dmean'
        torank            = abs(mean(data(classes==0,:), 'omitnan') - ...
                                mean(data(classes==1,:), 'omitnan'));
        [~, rankedList]   = sort(torank',    'descend'); %#ok<UDIM>
        [~, ranksOrdered] = sort(rankedList, 'ascend');
        
    case 'dmedian'
        torank            = abs(median(data(classes==0,:), 'omitnan') - ...
                                median(data(classes==1,:), 'omitnan'));
        [~, rankedList]   = sort(torank',    'descend'); %#ok<UDIM>
        [~, ranksOrdered] = sort(rankedList, 'ascend');
        
    case 'dstd'
        torank            = abs(std(data(classes==0,:), [], 'omitnan') - ...
                                std(data(classes==1,:), [], 'omitnan'));
        [~, rankedList]   = sort(torank',    'descend'); %#ok<UDIM>
        [~, ranksOrdered] = sort(rankedList, 'ascend');
        
    case 'correlation'
        torank            = corr(data, classes);
        [~, rankedList]   = sort(abs(torank),  'descend');
        [~, ranksOrdered] = sort(rankedList,   'ascend');
        
    case 'cohensd'
        term1             = mean(data(classes==0,:)) - mean(data(classes==1,:));
        term2             = ((sum(classes==0)-1).*std(data(classes==0,:))  +  ... 
                             (sum(classes==1)-1).*std(data(classes==1,:))) ./ ...
                             (length(classes) + 2);
        torank            = term1./term2;
        [~, rankedList]   = sort(abs(torank), 'descend');
        [~, ranksOrdered] = sort(rankedList,  'ascend');
        
    case 'hedgesg'
        term1             = mean(data(classes==0,:)) - mean(data(classes==1,:));
        term2             = ((sum(classes==0)-1).*std(data(classes==0,:)) +   ... 
                             (sum(classes==1)-1).*std(data(classes==1,:))) ./ ...
                             (length(classes) + 2);
        factor            = 1 - (3/(4*length(classes)-9));
        torank            = factor .* term1./term2;
        [~, rankedList]   = sort(abs(torank), 'descend');
        [~, ranksOrdered] = sort(rankedList,  'ascend');
        
    case 'glassd'
        term1             = mean(data(classes==0,:)) - mean(data(classes==1,:));
        term2             = std(data(classes==0,:));
        torank            = term1./term2;
        [~, rankedList]   = sort(abs(torank), 'descend');
        [~, ranksOrdered] = sort(rankedList,  'ascend');
        
    case 'silhouetteman'        
        % Number of observations in each class
        num_0 = sum(classes==0);
        num_1 = sum(classes==1);
        
        % Subset data
        sub_0 = data(classes==0,:);
        sub_1 = data(classes==1,:);
        
        % Initialize
        term_a0 = zeros(num_0, num_features);
        term_a1 = zeros(num_1, num_features);
        term_b0 = zeros(num_0, num_features);
        term_b1 = zeros(num_1, num_features);
        
        % Term a: average pairwise distance between all points within the
        % cluster, where average is taken with n-1
        % Term b: average pairwise distance between a point within the
        % cluster with all points in another cluster; note that min
        % operation is not needed as we assume that number of classes or
        % clusters is equal to two
        for pts = 1:num_0
            term_a0(pts,:) = sum(abs(sub_0(pts,:)  - sub_0))/(num_0 - 1);
            term_b0(pts,:) = mean(abs(sub_0(pts,:) - sub_1));
        end
        for pts = 1:num_1
            term_a1(pts,:) = sum(abs(sub_1(pts,:)    - sub_1))/(num_1 - 1);
            term_b1(pts,:) = mean(abs(sub_1(pts,:) - sub_0));
        end
        
        % Calculate Silhouette value for each point
        sil_0 = (term_b0 - term_a0)./max(term_a0, term_b0);
        sil_1 = (term_b1 - term_a1)./max(term_a1, term_b1);
        
        % Put values together and take the average
        torank = mean([sil_0; sil_1])';
        
        % Sort in descending order
        [~, rankedList]   = sort(abs(torank), 'descend');
        [~, ranksOrdered] = sort(rankedList,  'ascend');
        
    case 'silhouetteeuc'        
        % Number of observations in each class
        num_0 = sum(classes==0);
        num_1 = sum(classes==1);
        
        % Subset data
        sub_0 = data(classes==0,:);
        sub_1 = data(classes==1,:);
        
        % Initialize
        term_a0 = zeros(num_0, num_features);
        term_a1 = zeros(num_1, num_features);
        term_b0 = zeros(num_0, num_features);
        term_b1 = zeros(num_1, num_features);
        
        % Term a: average pairwise distance between all points within the
        % cluster, where average is taken with n-1
        % Term b: average pairwise distance between a point within the
        % cluster with all points in another cluster; note that min
        % operation is not needed as we assume that number of classes or
        % clusters is equal to two
        for pts = 1:num_0
            term_a0(pts,:) = sum(sqrt((sub_0(pts,:)  - sub_0).^2))/(num_0 - 1);
            term_b0(pts,:) = mean(sqrt((sub_0(pts,:) - sub_1).^2));
        end
        for pts = 1:num_1
            term_a1(pts,:) = sum(sqrt((sub_1(pts,:)  - sub_1).^2))/(num_1 - 1);
            term_b1(pts,:) = mean(sqrt((sub_1(pts,:) - sub_0).^2));
        end
        
        % Calculate Silhouette value for each point
        sil_0 = (term_b0 - term_a0)./max(term_a0, term_b0);
        sil_1 = (term_b1 - term_a1)./max(term_a1, term_b1);
        
        % Put values together and take the average
        torank = mean([sil_0; sil_1])';
        
        % Sort in descending order
        [~, rankedList]   = sort(abs(torank), 'descend');
        [~, ranksOrdered] = sort(rankedList,  'ascend');                
                
    otherwise
        error(['Unknown ranking method specified: ', rank_method]);
end