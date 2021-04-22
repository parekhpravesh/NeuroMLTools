function [sc_feat_matrix, sc_params] = feature_scaling(feat_matrix, method, sc_params)
% Function to scale a set of features, given a scaling method 
%% Inputs:
% feat_matrix:      [n x p] matrix with n samples and p features
% method:           should be one of the following (see Notes):
%                       * rescale
%                       * mean
%                       * std
% sc_params:        structure having scaling parameters to be applied to 
%                   feat_matrix depending on method (optional; see Notes):
%                       * if method == 'rescale', then should contain: 
%                               - min_val (minimum value per feature)
%                               - max_val (maximum value per feature)
%                       * if method == 'mean', then should contain: 
%                               - min_val  (minimum value per feature)
%                               - max_val  (maximum value per feature)
%                               - mean_val (mean value per feature)
%                       * if method == 'std', then should contain:
%                               - mean_val (mean value per feature)
%                               - std_val  (standard deviation per feature)
%                       
%% Outputs:
% sc_feat_matrix:   matrix of scaled features
% sc_params:        structure containing values used during feature scaling
% 
%% Defaults:
% method:           'rescale'
% sc_params:        [] (i.e. parameters are estimated)
% 
%% Notes:
% Method for feature scaling is based on the Wikipedia article on feature
% scaling: https://en.wikipedia.org/wiki/Feature_scaling
% 
% For each feature p in feat_matrix:
% * rescale
%       normalized_value = (x - min(x))/(max(x) - min(x))
% 
% * mean
%       normalized_value = (x - mean(x))/(max(x) - min(x))
% 
% * std
%       normalized_valie = (x - mean(x))/std(x)
% 
% If no scaling parameters are passed, the minimum, maximum, mean, and 
% standard deviation values (as needed) are calculated from the data; 
% these values are returned as a structure in sc_params; this is typically 
% the case when the input is the training data
% 
% If parameters are passed, these parameters are simply applied to the
% data; this is typically the case when applying the same feature scaling
% as training data to a test data; the same values are returned as
% sc_params
% 
% Variables in the sc_params structure should have the same number of 
% values as the number of features and should be ordered the same as the
% actual feature vector; sanity check is restricted to checking the number
% of features
% 
% If sc_params is passed as input, the output sc_params is the same as input
% 
%% Author(s):
% Parekh, Pravesh
% May 24, 2018
% April 22, 2021
% MBIAL

%% Validate inputs and assign defaults
% Check feature_matrix
if ~exist('feat_matrix', 'var') || isempty(feat_matrix)
    error('Feature matrix should be provided');
else
    % Number of features is the number of columns
    num_features = size(feat_matrix, 2);
end

% Check method
if ~exist('method', 'var') || isempty(method)
    method = 'rescale';
else
    method = lower(method);
end

% Check scaling_parameters
if ~exist('sc_params', 'var') || isempty(sc_params)
    to_scale = true;
else
    to_scale = false;
    % Check if required fields are present
    if strcmpi(method, 'rescale')
        if ~isfield(sc_params, 'min_val') || ...
           ~isfield(sc_params, 'max_val')
                error('min_val and max_val values should be provided');
        else
            % Check for sizes of scaling parameters
            if numel(sc_params.min_val) ~= num_features || ...
               numel(sc_params.max_val) ~= num_features
                error(['Mismatch between number of features: ', ...
                       num2str(num_features), ' and the scaling parameters']);
            end
        end
    else
        if strcmpi(method, 'mean')
            if ~isfield(sc_params, 'min_val') || ...
               ~isfield(sc_params, 'max_val') || ...
               ~isfield(sc_params, 'mean_val')
                    error('min_val, max_val, and mean_val should be provided');
            else
                % Check for sizes of scaling parameters
                if numel(sc_params.min_val)  ~= num_features || ...
                   numel(sc_params.max_val)  ~= num_features || ...
                   numel(sc_params.mean_val) ~= num_features
                        error(['Mismatch between number of features:' ,  ...
                               num2str(num_features), ' and the scaling parameters']);
                end
            end
        else
            if ~isfield(sc_params, 'mean_val') || ...
               ~isfield(sc_params, 'std_val')
                    error('mean_val and std_val should be provided');
            else
                % Check for sizes of scaling parameters
                if numel(sc_params.mean_val) ~= num_features || ...
                   numel(sc_params.std_val)  ~= num_features
                        error(['Mismatch between number of features ', ...
                               num2str(num_features), ' and the scaling parameters']);
                end
            end
        end
    end
end

%% Rescale each feature independently
switch(method)
    case 'rescale'
        if to_scale
            % Get minimum and maximum values
            sc_params.min_val = min(feat_matrix, [], 'omitnan');
            sc_params.max_val = max(feat_matrix, [], 'omitnan');
        end
        
        % Apply scaling
        sc_feat_matrix = bsxfun(@rdivide, bsxfun(@minus, feat_matrix,              sc_params.min_val), ...
                                                 bsxfun(@minus, sc_params.max_val, sc_params.min_val));

    case 'mean'
        if to_scale
            % Get minimum, maximum, and mean values
            sc_params.min_val  = min(feat_matrix,  [], 'omitnan');
            sc_params.max_val  = max(feat_matrix,  [], 'omitnan');
            sc_params.mean_val = mean(feat_matrix,     'omitnan');
        end
        
        % Apply scaling
        sc_feat_matrix = bsxfun(@rdivide, bsxfun(@minus, feat_matrix,              sc_params.mean_val), ...
                                                 bsxfun(@minus, sc_params.max_val, sc_params.min_val));

    case 'std'
        if to_scale
            % Get mean and standard deviation
            sc_params.mean_val = mean(feat_matrix,    'omitnan');
            sc_params.std_val  = std(feat_matrix, [], 'omitnan');
        end
        
        % Apply scaling
        sc_feat_matrix = bsxfun(@rdivide, bsxfun(@minus, feat_matrix, sc_params.mean_val), ...
                                                 sc_params.std_val);
                                             
    otherwise
        error(['Unknown feature scaling method provided: ', method]);
end