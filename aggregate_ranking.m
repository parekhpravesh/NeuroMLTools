function [ranks, aggregate] = aggregate_ranking(rank_mat, aggregation_method, rank_tie_breaker)
% Function to perform different types of rank aggregation, given a set of 
% feature ranking obtained by different ranking procedures
%% Inputs:
% rank_mat:             matrix of ranks where rows indicate the rank of
%                       that particular feature while columns indicate 
%                       ranking method; for example row=10, col=2 should 
%                       correspond to the rank of feature 10 obtained by
%                       ranking method 2 
% aggregation_method:   should be one of the following:
%                       *   'min' 
%                       *   'max'
%                       *   'mean'
%                       *   'median'
%                       *   'geomean'
%                       *   'roundrobin'
%                       *   'borda'
%                       *   'meanborda'
%                       *   'medianborda'
%                       *   'condorcet'
%                       *   'reciprocal'
% rank_tie_breaker:     true or false indicating if ties should be broken
%                       randomly or not
% 
%% Outputs:
% ranks:                a vector of feature numbers as per their ranking
%                       i.e., if first entry is 50, that means that feature
%                       number 50 was ranked the highest after median
%                       aggregate ranking
% aggregate:            a vector of aggregate ranks before any sorting; for
%                       example, the entry 50 would just be the computed
%                       aggregate ranking of feature 50; if aggregation
%                       method is Round Robin or condorcet, the aggregate 
%                       is empty
% 
%% Notes:
% Assume that that all features in the original feature matrix were ranked
% by all ranking methods
% 
% Aggregation method = 'min'
% --------------------------
%   * Divide each rank (within each column) by the total number of features
%   * For each feature, find the minimum rank across all ranking methods
%   * Sort the minimum rank in ascending order and assign this rank to that 
%     particular feature
%   * Break ties
% 
% Aggregation method = 'max'
% --------------------------
%   * Divide each rank (within each column) by the total number of features
%   * For each feature, find the maximum rank across all ranking methods
%   * Sort the maximum rank in ascending order and assign this rank to that 
%     particular feature
%   * Break ties
% 
% Aggregation method = 'mean'
% ---------------------------
%   * Divide each rank (within each column) by the total number of features
%   * For each feature, find the average rank across ranking methods
%   * Sort the mean in ascending order and assign this rank to that 
%     particular feature
%   * Break ties
% 
% Aggregation method = 'median'
% -----------------------------
%   * Divide each rank (within each column) by the total number of features
%   * For each feature, find the median rank across ranking methods
%   * Sort the median in ascending order and assign this rank to that 
%     particular feature
%   * Break ties
% 
% Aggregation method = 'geomean' (geometric mean)
% -----------------------------------------------
%   * Divide each rank (within each column) by the total number of features
%   * For each feature, find the product across ranking methods and then
%     take the nth root, where n = number of ranking methods
%   * Sort the geometric mean in ascending order and assign this rank to
%     that particular feature
%   * Break ties
%  
% Aggregation method = 'roundrobin' (Round Robin)
% -----------------------------------------------
%   * Sort each ranking column in ascending order such that the top
%     features in each ranking method are towards the top
%   * Assign a random ordering to the ranking methods
%   * First   feature of final list = first feature from the first list
%   * Second  feature of final list = first feature from the second list
%   * nth     feature of final list = first feature from the nth list
%   * (n+1)th feature of final list = second feature from the first list
%   * Proceed in this manner till all features are found
%   * Disregard any feature which was previously found
% 
% Aggregation method = 'borda'
% ----------------------------
%   * Assign a weight to each feature rank, for each ranking method, such
%     that:
%       - top feature:         weight = number of ranks
%       - second best feature: weight = number of ranks - 1
%       - worst rated feature: weight = 1
%   * Generate sum of weights across all ranking methods
%   * Sort by summed up weights in descending order
%   * Break ties randomly
% 
% Aggregation method = 'meanborda' (mean Borda)
% ---------------------------------------------
%   * Assign a weight to each feature rank, for each ranking method, such
%     that:
%       - top feature:         weight = number of ranks
%       - second best feature: weight = number of ranks - 1
%       - worst rated feature: weight = 1
%   * Generate average of weights across all ranking methods
%   * Sort by averaged weights in descending order
%   * Break ties randomly
% 
% Aggregation method = 'medianborda' (median Borda)
% -------------------------------------------------
%   * Assign a weight to each feature rank, for each ranking method, such
%     that:
%       - top feature:         weight = number of ranks
%       - second best feature: weight = number of ranks - 1
%       - worst rated feature: weight = 1
%   * Generate median of weights across all ranking methods
%   * Sort by median weights in descending order
%   * Break ties randomly
% 
% Aggregation method = 'condorcet' (Condorcet method)
% ---------------------------------------------------
%   * Create a preference matrix where each entry indicates the number of
%     times the given feature was ranked lower than all the other features
%     across ranking methods
%   * Calculate Copeland score by pairwise comparison of preference matrix:
%     a feature which dominates over the other, gets a value of 1 and a 
%     feature which has same preference as the other, gets a value of 0.5
%   * Calculate the Copeland result by summing the Copeland scores 
%   * Calculate the total number of times (across ranking methods) that a
%     particular feature was preferred over the other features
%   * Sort Copeland results in descending order to get feature ordering
%   * For ties in Copeland results, sort in descending order by total
%     number of times the feature was preferred
%   * If ties persist, break them randomly
% 
% Aggregation method = 'reciprocal'
% ---------------------------------
%   * Calculate the reciprocal fusion rank of a given feature as the sum of
%     reciprocals of (k + individual ranks) across ranking methods
%   * Sort in descending order of reciprocal fusion rank
%   * A value of k = 60 is used, as mentioned in Cormack et al (2009)
%   * Break ties randomly

% If rank_tie_breaker is false, then duplicate ranks are in the ascending 
% order of feature index; however, it is possible that more than one 
% feature gets the same rank; if rank_tie_breaker is true, then rank ties 
% are broken randomly
% 
% rank_tie_breaker is not applicable for Round Robin or Condorcet method
% 
% The description of Round Robin is based on the description in Dittman et
% al (2013); I assume that for Round Robin, it would be more meaningful to
% have the top ranking features towards the top of the list; otherwise, the
% results would not be top n features
% 
% The Condorcet method above is a simplistic implementation; it is possible
% to use graph theory methods like Schulze method to get a list of feature
% preferences
% 
%% Defaults:
% aggregation_method:   median
% rank_tie_breaker:     true
% 
%% References:
% Dittman D.J., Khoshgoftaar T.M., Wald R., and Napolitano A. 2013.
% Classification Performance of Rank Aggregation Techniques for 
% Ensemble Gene Selection. Proceedings of the Twenty-Sixth International 
% Florida Artificial Intelligence Research Society Conference
% https://www.aaai.org/ocs/index.php/FLAIRS/FLAIRS13/paper/viewFile/5893/6111
% 
% van Erp M., and Schomaker L. 2000. Variants of the Borda count method for
% combining ranked classifier hypotheses. In: L.R.B. Schomaker and L.G. 
% Vuurpijl (Eds.), Proceedings of the Seventh International Workshop on 
% Frontiers in Handwriting Recognition,September 11-13 2000, Amsterdam, 
% ISBN 90-76942-01-3, Nijmegen: International Unipen Foundation, pp 443-452
% http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.30.1166
% 
% Cormack G.V., Clarke, C.L.A., Büttcher S. 2009. Reciprocal Rank Fusion
% outperforms Condorcet and individual Rank Learning Methods. 
% SIGIR '09: Proceedings of the 32nd international ACM SIGIR conference on 
% Research and development in information retrieval Pages 758–759
% https://doi.org/10.1145/1571941.1572114
% http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.150.2291&rep=rep1&type=pdf
% 
% % Find duplicate rankings using Guillaume's solution found at:
% https://in.mathworks.com/matlabcentral/answers/388695
% 
%% See also:
% Nejc Ilc (2021). Rank aggregation 
% (https://www.mathworks.com/matlabcentral/fileexchange/41835-rank-aggregation), 
% MATLAB Central File Exchange. Retrieved April 17, 2021. 
% 
%% Author(s)
% Parekh, Pravesh
% April 17, 2021
% MBIAL

%% Check inputs
% Check rank_mat
if ~exist('rank_mat', 'var') || isempty(rank_mat)
    error('Please provide a ranking matrix');
else
    [num_ranks, num_methods] = size(rank_mat);
end

% Check aggregation_method
if ~exist('aggregation_method', 'var') || isempty(aggregation_method)
    aggregation_method = 'median';
end

% Check rank_tie_breaker
if ~exist('rank_tie_breaker', 'var') || isempty(rank_tie_breaker)
    rank_tie_breaker = true;
else
     if ~islogical(rank_tie_breaker)
         error('rank_tie_breaker should be either true or false');
     end
end

%% Aggregate ranks
switch aggregation_method
    case 'min'        
        % Normalize ranks, find minimum across ranking methods, sort
        norm_ranks  = rank_mat./num_ranks;
        aggregate   = min(norm_ranks, [], 2);
        [~, ranks]  = sort(aggregate, 'ascend');
        
    case 'max'
        % Normalize ranks, find maximum across ranking methods, sort
        norm_ranks  = rank_mat./num_ranks;
        aggregate   = max(norm_ranks, [], 2);
        [~, ranks]  = sort(aggregate, 'ascend');

    case 'mean'
        % Normalize ranks, find mean across ranking methods, sort
        norm_ranks  = rank_mat./num_ranks;
        aggregate   = mean(norm_ranks, 2);
        [~, ranks]  = sort(aggregate, 'ascend');
        
    case 'median'
        % Normalize ranks, find median across ranking methods, sort
        norm_ranks  = rank_mat./num_ranks;
        aggregate   = median(norm_ranks, 2);
        [~, ranks]  = sort(aggregate, 'ascend');
        
    case 'geomean'
        % Normalize ranks, find geometric mean across ranking methods, sort
        norm_ranks  = rank_mat./num_ranks;
        aggregate   = nthroot(prod(norm_ranks, 2), num_methods);
        [~, ranks]  = sort(aggregate, 'ascend');
        
    case 'roundrobin'
        % Sort each column in the rank_mat
        for fs = 1:num_methods
            [~, rank_mat(:,fs)] = sort(rank_mat(:,fs), 'ascend');
        end
        
        % Randomly reorder the rank_mat columns
        rank_mat = rank_mat(:,randperm(4,4));
        
        % Initialize
        ranks     = [];
        aggregate = [];
        loc       = 1;
        
        % Loop over which value to pick
        for pick = 1:num_ranks
            
            % Loop over lists
            for list = 1:num_methods
                
                % Pick a value
                tmp = rank_mat(pick, list);
                
                % Check if value is already part of ranks
                if sum(tmp == ranks)
                    % Do nothing
                else
                    % Add feature to ranks
                    ranks(loc, 1) = tmp; %#ok<AGROW>
                    loc           = loc + 1;
                end
            end
        end
        
    case 'borda'
        % Create weights and add them together
        tmp_weights = ((num_ranks:-1:1)');
        weights     = tmp_weights(rank_mat);
        aggregate   = sum(weights,2);
        
        % Sort in descending order
        [~, ranks] = sort(aggregate, 'descend');
        
    case 'meanborda'
        % Create weights and take their average
        tmp_weights = ((num_ranks:-1:1)');
        weights     = tmp_weights(rank_mat);
        aggregate   = mean(weights,2);
        
        % Sort in descending order
        [~, ranks] = sort(aggregate, 'descend');

    case 'medianborda'
        % Create weights and take their median
        tmp_weights = ((num_ranks:-1:1)');
        weights     = tmp_weights(rank_mat);
        aggregate   = median(weights,2);
        
        % Sort in descending order
        [~, ranks] = sort(aggregate, 'descend');
        
    case 'condorcet'
        % Create a preference matrix where each row is the number of times
        % a given feature was ranked better than remaining features
        pref_mat = zeros(num_ranks, num_ranks);
        for fs = 1:num_ranks
            pref_mat(fs,:) = sum(rank_mat(fs,:) < rank_mat, 2)';
        end
        
        % Create a Copeland result matrix
        copeland = zeros(num_ranks, num_ranks);
        for fs1 = 1:num_ranks
            for fs2 = 1:num_ranks
                if pref_mat(fs1,fs2) > pref_mat(fs2,fs1)
                    copeland(fs1,fs2) = 1;
                else
                    if pref_mat(fs1,fs2) == pref_mat(fs2,fs1)
                        copeland(fs1,fs2) = 0.5;
                    end
                end
            end
        end
        % Convert diagonals back to zero
        copeland = copeland - diag(diag(copeland));
        
        % Find overall Copeland result vector
        aggregate = sum(copeland,2);
        
        % Find the total number of rankers which preferred a feature over
        % other features
        tot_pref = sum(pref_mat,2);
        
        % Sort in descending order based on Copeland results
        [~, ranks] = sort(aggregate, 'descend');
        
        % Rearrange aggreagate and tot_pref in the same order as ranks
        tot_pref  = tot_pref(ranks);
        aggregate = aggregate(ranks);
        
        % Find unique Copeland scores
        uq_Copeland = unique(aggregate);
        
        % Go over every Copeland score and handle duplicates
        for uq = 1:length(uq_Copeland)
            tmp_count = aggregate == uq_Copeland(uq);
            if sum(tmp_count) > 1
                
                % Duplicate found; get their total number of rankers
                loc_dup = find(aggregate == uq_Copeland(uq));
                tmp_tot = tot_pref(loc_dup);
                
                % Sort in descending order and rearrange ranks and tot_pref
                [~, b]            = sort(tmp_tot, 'descend');
                ranks(loc_dup)    = ranks(loc_dup(b));
                tot_pref(loc_dup) = tot_pref(loc_dup(b));
                
                % Check if there were duplicates in tmp_tot
                tmp_uq  = unique(tmp_tot);
                for tt  = 1:length(tmp_uq)
                    tmp = tmp_tot == tmp_uq(tt);
                    if sum(tmp) > 1
                        loc_dDup = find(tmp_tot == tmp_uq(tt));
                        tmp_ord  = randperm(length(loc_dDup), length(loc_dDup));
                        ranks(loc_dup(loc_dDup))    = ranks(loc_dup(loc_dDup(tmp_ord)));
                        tot_pref(loc_dup(loc_dDup)) = tot_pref(loc_dup(loc_dDup(tmp_ord)));
                    end
                end
            end
        end
        aggregate = [];
        
    case 'reciprocal'
        aggregate  = sum(1./(60 + rank_mat),2);
        [~, ranks] = sort(aggregate, 'descend');
        
    otherwise
        error(['Unknown aggregation method specified: ', aggregation_method]);
end

%% Resolve rank ties, if needed
% Find duplicate rankings using Guillaume's solution
% https://in.mathworks.com/matlabcentral/answers/388695
if rank_tie_breaker
    [uvalues, ~, uid]       = unique(aggregate(:));
    count                   = accumarray(uid, 1);                                        %easiest way to calculate the histogram of uvalues
    linindices              = accumarray(uid, (1:numel(aggregate))', [], @(idx) {idx});  %split linear indices according to uid
    valwhere                = [num2cell(uvalues), linindices];                           %concatenate
    valwhere(count == 1, :) = [] ;                                                       %remove count of 1
    
    % Go over each duplicate and resolve them randomly
    for rep = 1:size(valwhere,1)
        % Find these locations in rank
        loc = [];
        for tmp = 1:length(valwhere{rep,2})
            loc(tmp,1) = find(ranks == valwhere{rep,2}(tmp));  %#ok<AGROW>
        end
        % Permute
        ranks(loc) = ranks(loc(randperm(length(loc), length(loc))));
    end
end