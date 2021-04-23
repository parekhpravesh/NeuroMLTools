function globalSurf = get_cat_globalsurf(surfaces, dir_output, toWrite)
% Function to compile whole brain surface summary measures from CAT 
% segmentation
%% Inputs:
% surfaces:         full path to surf directory containing various surface
%                   data files (such as cortical thickness) OR cell type
%                   with a list of full paths to LH PBT files
% dir_output:       full path to where results should be saved
% toWrite:          true or false indicating if the csv file should be
%                   written out
% 
%% Output:
% globalSurf is a table type variable that contains the following columns:
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
% If toWrite is true, a csv file named 'SurfMeasures_WholeBrain' is written 
% in dir_output that contains the above columns
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
% dir_output:   same as first entry in surfaces or pwd
% toWrite:      true
% 
%% Author(s):
% Parekh, Pravesh
% April 21, 2021
% MBIAL

%% Check inputs
% Check surfaces
if ~exist('surfaces', 'var') || isempty(surfaces)
    error('Please provide full path to the surf directory or a list of PBT files to work on');
else
    if iscell(surfaces)
        num_surfaces = length(surfaces);
        % Make sure all PBT files exist
        f = @(x) exist(x, 'file');
        if sum(logical(cellfun(f, surfaces))) ~= num_surfaces
            error('One or more surface files not found');
        end
    else
        if ~exist(surfaces, 'dir')
            error(['Unable to find: ', surfaces]);
        else
            tmp_surfaces = dir(fullfile(surfaces, 'lh.pbt.*'));
            tmp_surfaces = {tmp_surfaces(:).name}';
            surfaces     = fullfile(surfaces, tmp_surfaces);
            num_surfaces = length(surfaces);
        end
    end
end

% Check dir_output
if ~exist('dir_output', 'var') || isempty(dir_output)
    % Check if full path is present in the first file
    if isempty(fileparts(surfaces{1}))
        outname = fullfile(pwd, 'SurfMeasures_WholeBrain.csv');
    else
        outname = fullfile(fileparts(surfaces{1}), 'SurfMeasures_WholeBrain.csv');
    end
else
    if ~exist(dir_output, 'dir')
        mkdir(dir_output);
    end
    outname = fullfile(dir_output, 'SurfMeasures_WholeBrain.csv');
end

% Check toWrite
if ~exist('toWrite', 'var') || isempty(toWrite)
    toWrite = true;
else
    if ~islogical(toWrite)
        error('toWrite should be either true or false');
    end
end

%% Compile measures
globalSurf = cell(num_surfaces, 37);

% Loop over each subject file and summarize values
for subjs  = 1:num_surfaces
        
    % Get subject ID
    [~, tmp_name]       = fileparts(regexprep(strrep(surfaces{subjs}, 'lh.pbt.', ''), {'_T1w_RAS', '_T1w'}, ''));
    globalSurf{subjs,1} = tmp_name;
    
    % Get PBT info
    dat_LH                  = cat_io_FreeSurfer('read_surf_data', surfaces{subjs});
    dat_RH                  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'rh.pbt.'));
    globalSurf{subjs, 2}    = mean([dat_LH; dat_RH]);
    globalSurf{subjs, 3}    = std([dat_LH;  dat_RH]);
    globalSurf{subjs, 14}   = mean(dat_LH);
    globalSurf{subjs, 15}   = std(dat_LH);
    globalSurf{subjs, 16}   = mean(dat_RH);
    globalSurf{subjs, 17}   = std(dat_RH);
    
    % Get cortical thickness
    try
        dat_LH                  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'lh.thickness.'));
        dat_RH                  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'rh.thickness.'));
        globalSurf{subjs, 4}    = mean([dat_LH; dat_RH]);
        globalSurf{subjs, 5}    = std([dat_LH;  dat_RH]);
        globalSurf{subjs, 18}   = mean(dat_LH);
        globalSurf{subjs, 19}   = std(dat_LH);
        globalSurf{subjs, 20}   = mean(dat_RH);
        globalSurf{subjs, 21}   = std(dat_RH);
    catch
        warning(['Cortical thickness file not found for: ', tmp_name]);
    end
    
    % Get gyrification
    try
        dat_LH                  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'lh.gyrification.'));
        dat_RH                  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'rh.gyrification.'));
        globalSurf{subjs, 6}    = mean([dat_LH; dat_RH]);
        globalSurf{subjs, 7}    = std([dat_LH;  dat_RH]);
        globalSurf{subjs, 22}   = mean(dat_LH);
        globalSurf{subjs, 23}   = std(dat_LH);
        globalSurf{subjs, 24}   = mean(dat_RH);
        globalSurf{subjs, 25}   = std(dat_RH);
    catch
        warning(['Gyrification file not found for: ', tmp_name]);
    end
    
    % Get sulcal depth
    try
        dat_LH                  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'lh.depth.'));
        dat_RH                  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'rh.depth.'));
        globalSurf{subjs, 8}    = mean([dat_LH; dat_RH]);
        globalSurf{subjs, 9}    = std([dat_LH;  dat_RH]);
        globalSurf{subjs, 26}   = mean(dat_LH);
        globalSurf{subjs, 27}   = std(dat_LH);
        globalSurf{subjs, 28}   = mean(dat_RH);
        globalSurf{subjs, 29}   = std(dat_RH);
    catch
        warning(['Depth file not found for: ', tmp_name]);
    end
    
    % Get fractal dimension
    try
        dat_LH                  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'lh.fractaldimension.'));
        dat_RH                  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'rh.fractaldimension.'));
        globalSurf{subjs, 10}   = mean([dat_LH; dat_RH]);
        globalSurf{subjs, 11}   = std([dat_LH;  dat_RH]);
        globalSurf{subjs, 30}   = mean(dat_LH);
        globalSurf{subjs, 31}   = std(dat_LH);
        globalSurf{subjs, 32}   = mean(dat_RH);
        globalSurf{subjs, 33}   = std(dat_RH);
    catch
        warning(['Depth file not found for: ', tmp_name]);
    end    
    
    % Get Toro's gyrification index
    try
        tmp_LH                  = dir(fullfile(fileparts(surfaces{subjs}), ['lh.toroGI*.', tmp_name, '*']));
        tmp_RH                  = dir(fullfile(fileparts(surfaces{subjs}), ['rh.toroGI*.', tmp_name, '*']));
        dat_LH                  = cat_io_FreeSurfer('read_surf_data', fullfile(fileparts(surfaces{subjs}), tmp_LH(1).name));
        dat_RH                  = cat_io_FreeSurfer('read_surf_data', fullfile(fileparts(surfaces{subjs}), tmp_RH(1).name));
        globalSurf{subjs, 12}   = mean([dat_LH; dat_RH]);
        globalSurf{subjs, 13}   = std([dat_LH;  dat_RH]);
        globalSurf{subjs, 34}   = mean(dat_LH);
        globalSurf{subjs, 35}   = std(dat_LH);
        globalSurf{subjs, 36}   = mean(dat_RH);
        globalSurf{subjs, 37}   = std(dat_RH);
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
globalSurf = cell2table(globalSurf, 'VariableNames', var_names);

% Write table, if required
if toWrite
    writetable(globalSurf, outname);
end