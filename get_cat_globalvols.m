function get_cat_globalvols(dir_reports, dir_output)
% Function to compile whole brain volumetric summary measures from CAT 
% segmentation
%% Inputs:
% dir_reports:      full path to report directory containing cat_*.mat files
% dir_output:       full path to where results should be saved
% 
%% Output:
% An csv file named 'Volumes_WholeBrain' is written in dir_output that
% contains the following information:
%   * SubjectID
%   * Total GM volume
%   * Total WM volume
%   * Total CSF volume
%   * Total WMH volume
%   * Total GM percentage relative to TIV
%   * Total WM percentage relative to TIV
%   * Total CSF percentage relative to TIV
%   * Total WMH percentage relative to TIV
%   * Total intracranial volume (TIV)
% 
%% Defaults:
% dir_output:   same as dir_reports
% 
%% Author(s):
% Parekh, Pravesh
% April 21, 2021
% MBIAL

%% Check inputs
% Check dir_reports
if ~exist('dir_reports', 'var') || isempty(dir_reports)
    error('Please provide full path to the reports directory');
else
    if ~exist(dir_reports, 'dir')
        error(['Unable to find: ', dir_reports]);
    end
end

% Check dir_output
if ~exist('dir_output', 'var') || isempty(dir_output)
    dir_output = dir_reports;
else
    if ~exist(dir_output, 'dir')
        mkdir(dir_output);
    end
end

%% Compile results
list_mat = dir(fullfile(dir_reports, 'cat_*.mat'));
if isempty(list_mat)
    error(['Unable to find any cat_*.mat files in: ', dir_reports]);
else
    rel_vol  = zeros(length(list_mat),4);
    abs_vol  = zeros(length(list_mat),4);
    TIV      = zeros(length(list_mat),1);
    
    % Loop over each mat file and get values
    for mat  = 1:length(list_mat)
        load(fullfile(dir_reports, list_mat(mat).name), 'S');
        tmp              = S.subjectmeasures.vol_rel_CGW*100;
        rel_vol(mat,1:4) = tmp([2,3,1,4]);
        abs_vol(mat,1:4) = S.subjectmeasures.vol_abs_CGW([2,1,3,4]);
        TIV(mat,1)       = S.subjectmeasures.vol_TIV;
    end
    
    % Put together as a table
    res         = cell(length(list_mat),10);
    res(:,1)    = regexprep({list_mat(:).name}', {'cat_', '_T1w_RAS', '_T1w', '.mat'}, '');
    res(:,2:5)  = num2cell(abs_vol);
    res(:,6:9)  = num2cell(rel_vol);
    res(:,10)   = num2cell(TIV);
    res         = cell2table(res, 'VariableNames', {'SubjectID'; 'GM'; 'WM'; 'CSF'; 'WMH'; 'RelGM'; 'RelWM'; 'RelCSF'; 'RelWMH'; 'TIV'});
    writetable(res, fullfile(dir_output, 'Volumes_WholeBrain.csv'));
end