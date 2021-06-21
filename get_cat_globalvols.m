function globalVols = get_cat_globalvols(reports, dir_output, toWrite)
% Function to compile whole brain volumetric summary measures from CAT 
% segmentation
%% Inputs:
% reports:          full path to report directory containing cat_*.mat
%                   files OR a cell type having full paths to cat_*.mat
%                   files (rows are filenames)
% dir_output:       full path to where results should be saved
% toWrite:          true or false indicating if the csv file should be
%                   written out
% 
%% Outputs:
% globalVols is a table that contains the following columns:
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
% If toWrite is true, then a csv file named 'Volumes_WholeBrain' is written 
% in dir_output that contains the above columns
% 
%% Defaults:
% dir_output:       same as first entry in reports OR pwd
% toWrite:          true
% 
%% Author(s):
% Parekh, Pravesh
% April 21, 2021
% MBIAL

%% Check inputs
% Check reports
if ~exist('reports', 'var') || isempty(reports)
    error('Please provide full path to the reports directory or a list of reports to work on');
else
    if iscell(reports)
        num_reports = length(reports);
        % Make sure all report files exist
        f = @(x) exist(x, 'file');
        if sum(logical(cellfun(f, reports))) ~= num_reports
            error('One or more report files not found');
        end
    else
        if ~exist(reports, 'dir')
            error(['Unable to find: ', reports]);
        else
            tmp_reports = dir(fullfile(reports, 'cat_*.mat'));
            tmp_reports = {tmp_reports(:).name}';
            reports     = fullfile(reports, tmp_reports);
            num_reports = length(reports);
        end
    end
end

% Check dir_output
if ~exist('dir_output', 'var') || isempty(dir_output)
    % Check if full path is present in the first file
    if isempty(fileparts(reports{1}))
        outname = fullfile(pwd, 'Volumes_WholeBrain.csv');
    else
        outname = fullfile(fileparts(reports{1}), 'Volumes_WholeBrain.csv');
    end
else
    if ~exist(dir_output, 'dir')
        mkdir(dir_output);
    end
    outname = fullfile(dir_output, 'Volumes_WholeBrain.csv');
end

% Check toWrite
if ~exist('toWrite', 'var') || isempty(toWrite)
    toWrite = true;
else
    if ~islogical(toWrite)
        error('toWrite should be either true or false');
    end
end

%% Compile results
% Initialize
rel_vol    = zeros(num_reports,4);
abs_vol    = zeros(num_reports,4);
TIV        = zeros(num_reports,1);
header     = {'SubjectID'; 'GM'; 'WM'; 'CSF'; 'WMH'; 'RelGM'; 'RelWM'; 'RelCSF'; 'RelWMH'; 'TIV'};
globalVols = cell(num_reports, length(header));
locError   = [];

% Loop over each mat file and get values
for report = 1:num_reports
    
    % Load report
    load(reports{report}, 'S');
    
    % Get subject ID
    [~, globalVols{report,1}] = fileparts(regexprep(reports{report}, {'cat_', '_T1w_RAS', '_T1w', '.mat'}, ''));
    
    % Skip subject if there was an error
    if isfield(S, 'error')
        warning(['Error in report: ', reports{report}, '; skipping this report']);
        locError = [locError; report]; %#ok<AGROW>
        continue;
    end    

    tmp                 = S.subjectmeasures.vol_rel_CGW*100;
    rel_vol(report,1:4) = tmp([2,3,1,4]);
    abs_vol(report,1:4) = S.subjectmeasures.vol_abs_CGW([2,3,1,4]);
    TIV(report,1)       = S.subjectmeasures.vol_TIV;
end

% Put together as a table
globalVols(:,2:5)    = num2cell(abs_vol);
globalVols(:,6:9)    = num2cell(rel_vol);
globalVols(:,10)     = num2cell(TIV);
globalVols           = cell2table(globalVols, 'VariableNames', header);

% Remove rows that had errors
globalVols(locError,:) = [];

% Write table, if required
if toWrite
    writetable(globalVols, outname);
end