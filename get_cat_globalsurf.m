function get_cat_globalsurf(dir_surf, dir_output)
% Function to compile whole brain surface summary measures from CAT 
% segmentation
%% Inputs:
% dir_surf:         full path to surf directory containing various surface
%                   data files (such as cortical thickness)
% dir_output:       full path to where results should be saved
% 
%% Output:
% An csv file named 'SurfMeasures_WholeBrain' is written in dir_output that
% contains the following information:
%   * SubjectID
%   * Mean PBT cortical thickness
%   * Standard deviation of PBT cortical thickness 
% 
%   * Mean cortical thickness
%   * Standard deviation of cortical thickness
% 
%   * Mean gyrification
%   * Standard deviation of gyrification
% 
%   * Mean sulcal depth
%   * Standard deviation of sulcal depth
% 
%   * Mean fractal dimension
%   * Standard deviation of fractal dimension
% 
%   * Mean Toro's gyrification index
%   * Standard deviation of Toro's gyrification index
% 
%   * Mean PBT cortical thickness [left hemisphere]
%   * Standard deviation PBT cortical thickness [left hemisphere]
%   * Mean PBT cortical thickness [right hemisphere]
%   * Standard deviation PBT cortical thickness [right hemisphere]
% 
%   * Mean cortical thickness [left hemisphere]
%   * Standard deviation cortical thickness [left hemisphere]
%   * Mean cortical thickness [right hemisphere]
%   * Standard deviation cortical thickness [right hemisphere]
% 
%   * Mean gyrification [left hemisphere]
%   * Standard deviation gyrification [left hemisphere]
%   * Mean gyrification [right hemisphere]
%   * Standard deviation gyrification [right hemisphere]
% 
%   * Mean sulcal depth [left hemisphere]
%   * Standard deviation sulcal depth [left hemisphere]
%   * Mean sulcal depth [right hemisphere]
%   * Standard deviation sulcal depth [right hemisphere]
% 
%   * Mean fractal dimension [left hemisphere]
%   * Standard deviation fractal dimension [left hemisphere]
%   * Mean fractal dimension [right hemisphere]
%   * Standard deviation fractal dimension [right hemisphere]
% 
%   * Mean Toro's gyrification index [left hemisphere]
%   * Standard deviation Toro's gyrification index [left hemisphere]
%   * Mean Toro's gyrification index [right hemisphere]
%   * Standard deviation Toro's gyrification index [right hemisphere]
% 
%% Notes:
% If any of the surface measure files are missing, those measures are
% skipped from compilation
% 
% Assumes the following naming convention:
% Projection based cortical thickness:  [l/r]h.pbt.*
% Cortical thickness:                   [l/r]h.thickness.*
% Gyrification:                         [l/r]h.gyrification.*
% Sulcal depth:                         [l/r]h.depth.*
% Fractal dimension:                    [l/r]h.fractaldimension.*
% Toro's gyrification index:            [l/r]h.toroGI*.*
% 
% Note that CAT reports use PBT for reporting mean and standard deviation
% of cortical thickness
% 
%% Defaults:
% dir_output:   same as dir_surf
% 
%% Author(s):
% Parekh, Pravesh
% April 21, 2021
% MBIAL

%% Check inputs
% Check dir_reports
if ~exist('dir_surf', 'var') || isempty(dir_surf)
    error('Please provide full path to the surf directory');
else
    if ~exist(dir_surf, 'dir')
        error(['Unable to find: ', dir_surf]);
    end
end

% Check dir_output
if ~exist('dir_output', 'var') || isempty(dir_output)
    dir_output = dir_surf;
else
    if ~exist(dir_output, 'dir')
        mkdir(dir_output);
    end
end

%% Make a list of subjects 
list_pbt = dir(fullfile(dir_surf, 'lh.pbt.*'));
if isempty(list_pbt)
    error(['Unable to find any left hemisphere PBT files in: ', dir_surf, '; aborting']);
end

%% Compile measures
res = cell(length(list_pbt),37);

% Loop over each subject file and summarize values
for subjs  = 1:length(list_pbt)
        
    % Get subject ID
    tmp_name     = regexprep(strrep(list_pbt(subjs).name, 'lh.pbt.', ''), {'_T1w_RAS', '_T1w'}, '');
    res{subjs,1} = tmp_name;
    
    % Get PBT info
    dat_LH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, list_pbt(subjs).name));
    dat_RH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, strrep(list_pbt(subjs).name, 'lh.pbt.', 'rh.pbt.')));
    res{subjs, 2}   = mean([dat_LH; dat_RH]);
    res{subjs, 3}   = std([dat_LH;  dat_RH]);
    res{subjs, 14}  = mean(dat_LH);
    res{subjs, 15}  = std(dat_LH);
    res{subjs, 16}  = mean(dat_RH);
    res{subjs, 17}  = std(dat_RH);
    
    % Get cortical thickness
    try
        dat_LH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, strrep(list_pbt(subjs).name, 'lh.pbt.', 'lh.thickness.')));
        dat_RH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, strrep(list_pbt(subjs).name, 'lh.pbt.', 'rh.thickness.')));
        res{subjs, 4}   = mean([dat_LH; dat_RH]);
        res{subjs, 5}   = std([dat_LH;  dat_RH]);
        res{subjs, 18}  = mean(dat_LH);
        res{subjs, 19}  = std(dat_LH);
        res{subjs, 20}  = mean(dat_RH);
        res{subjs, 21}  = std(dat_RH);
    catch
        warning(['Cortical thickness file not found for: ', tmp_name]);
    end
    
    % Get gyrification
    try
        dat_LH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, strrep(list_pbt(subjs).name, 'lh.pbt.', 'lh.gyrification.')));
        dat_RH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, strrep(list_pbt(subjs).name, 'lh.pbt.', 'rh.gyrification.')));
        res{subjs, 6}   = mean([dat_LH; dat_RH]);
        res{subjs, 7}   = std([dat_LH;  dat_RH]);
        res{subjs, 22}  = mean(dat_LH);
        res{subjs, 23}  = std(dat_LH);
        res{subjs, 24}  = mean(dat_RH);
        res{subjs, 25}  = std(dat_RH);
    catch
        warning(['Gyrification file not found for: ', tmp_name]);
    end
    
    % Get sulcal depth
    try
        dat_LH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, strrep(list_pbt(subjs).name, 'lh.pbt.', 'lh.depth.')));
        dat_RH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, strrep(list_pbt(subjs).name, 'lh.pbt.', 'rh.depth.')));
        res{subjs, 8}   = mean([dat_LH; dat_RH]);
        res{subjs, 9}   = std([dat_LH;  dat_RH]);
        res{subjs, 26}  = mean(dat_LH);
        res{subjs, 27}  = std(dat_LH);
        res{subjs, 28}  = mean(dat_RH);
        res{subjs, 29}  = std(dat_RH);
    catch
        warning(['Depth file not found for: ', tmp_name]);
    end
    
    % Get fractal dimension
    try
        dat_LH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, strrep(list_pbt(subjs).name, 'lh.pbt.', 'lh.fractaldimension.')));
        dat_RH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, strrep(list_pbt(subjs).name, 'lh.pbt.', 'rh.fractaldimension.')));
        res{subjs, 10}  = mean([dat_LH; dat_RH]);
        res{subjs, 11}  = std([dat_LH;  dat_RH]);
        res{subjs, 30}  = mean(dat_LH);
        res{subjs, 31}  = std(dat_LH);
        res{subjs, 32}  = mean(dat_RH);
        res{subjs, 33}  = std(dat_RH);
    catch
        warning(['Depth file not found for: ', tmp_name]);
    end    
    
    % Get Toro's gyrification index
    try
        tmp_LH          = dir(fullfile(dir_surf, ['lh.toroGI*.', tmp_name, '*']));
        tmp_RH          = dir(fullfile(dir_surf, ['rh.toroGI*.', tmp_name, '*']));
        dat_LH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, tmp_LH(1).name));
        dat_RH          = cat_io_FreeSurfer('read_surf_data', fullfile(dir_surf, tmp_RH(1).name));
        res{subjs, 12}  = mean([dat_LH; dat_RH]);
        res{subjs, 13}  = std([dat_LH;  dat_RH]);
        res{subjs, 34}  = mean(dat_LH);
        res{subjs, 35}  = std(dat_LH);
        res{subjs, 36}  = mean(dat_RH);
        res{subjs, 37}  = std(dat_RH);
    catch
        warning(['Toro GI files not found for: ', tmp_name]);
    end        
end

% Put together as a table
var_names = {'SubjectID';                                               ...
             'Mean_PBT';                    'SD_PBT';                   ...
             'Mean_Thickness';              'SD_Thickness';             ...
             'Mean_Gyrification';           'SD_Gyrification';          ...
             'Mean_Depth';                  'SD_Depth';                 ...
             'Mean_FractalDimension';       'SD_FractalDimension';      ...
             'Mean_ToroGI';                 'SD_ToroGI';                ...
             'Mean_LH_PBT';                 'SD_LH_PBT';                ...
             'Mean_RH_PBT';                 'SD_RH_PBT';                ...
             'Mean_LH_Thickness';           'SD_LH_Thickness';          ...
             'Mean_RH_Thickness';           'SD_RH_Thickness';          ...
             'Mean_LH_Gyrification';        'SD_LH_Gyrification';       ...
             'Mean_RH_Gyrification';        'SD_RH_Gyrification';       ...
             'Mean_LH_Depth';               'SD_LH_Depth';              ...
             'Mean_RH_Depth';               'SD_RH_Depth';              ...
             'Mean_LH_FractalDimension';    'SD_LH_FractalDimension';   ...
             'Mean_RH_FractalDimension';    'SD_RH_FractalDimension';   ...
             'Mean_LH_ToroGI';              'SD_LH_ToroGI';             ...
             'Mean_RH_ToroGI';              'SD_RH_ToroGI'};
             
res = cell2table(res, 'VariableNames', var_names);
writetable(res, fullfile(dir_output, 'SurfMeasures_WholeBrain.csv'));
end