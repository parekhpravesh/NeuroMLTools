function [residuals, coeff] = regress_covariates(matrix, cov_mat, coeff)
% Function to regress covariates from each column of a matrix
%% Inputs:
% matrix:       [n x p] matrix where n are the number of observations and p 
%               are the number of columns; regression is performed on each 
%               column p
% cov_mat:      [n x k] matrix where n are the number of observations and k 
%               are the number of covariates
% coeff:        (optional) a set of estimated coefficients to use
% 
%% Outputs:
% residuals:    [n x p] matrix of residuals after regression
% coeff:        [1 x k] estimated linear regression coefficients
% 
%% Notes:
% For a given column p of matrix, we fit a linear regression model using
% the covariates and then calculate the residuals for that column p; the 
% coefficients for each covariate is also saved
% 
% If coeff is provided, then the coefficients are used to calculate
% residuals rather than fitting a linear regression model (useful for
% 'applying' a regression model to test data)
% 
% Automatically adds a constant term as the last column of cov_mat
% 
%% Author(s):
% Parekh, Pravesh
% January 23, 2020
% April 22, 2021
% MBIAL

%% Check and parse inputs
% Check feat_mat
if ~exist('matrix', 'var') || isempty(matrix)
    error('Please provide the matrix to perform regression on');
else
    [num_samples, num_feat] = size(matrix);
end

% Check cov_mat
if ~exist('cov_mat', 'var') || isempty(cov_mat)
    error('Please provide the covariates to regress out');
else
    [nsample, num_cov] = size(cov_mat);
    if num_samples ~= nsample
        error('Mismatch between number of observations in matrix and covariate matrix');
    end
    % Adding constant term
    cov_mat(:,end+1) = ones(nsample,1);
    num_cov          = num_cov + 1;
end

% Check coeff
if ~exist('coeff', 'var') || isempty(coeff)
    to_estimate = true;
else
    to_estimate = false;
    [nfeat, ncoeff] = size(coeff);
    if nfeat ~= num_feat
        error('Mismatch between number of columns in matrix and number of coefficients');
    end
    if ncoeff ~= num_cov
        error('Mismatch between number of covariates and number of coefficients');
    end
end

%% Initialize
if to_estimate
    coeff = zeros(num_feat, num_cov);
end
residuals = zeros(num_samples, num_feat);

%% Estimate coefficients, if necessary
if to_estimate
    for feat = 1:num_feat
        coeff(feat, 1:num_cov) = cov_mat\matrix(:,feat);
    end
end

%% Calculate residuals
for feat = 1:num_feat
    yhat = sum(bsxfun(@times, coeff(feat,1:num_cov), cov_mat),2);
    residuals(:,feat) = matrix(:,feat) - yhat;
end