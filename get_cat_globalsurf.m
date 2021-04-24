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

%% Handle Toro's gyrification indices and initialize
% Work out the number of Toro's gyrificiation indices that might be present
% using the first subject
tmp_Toro = dir(strrep(surfaces{1}, 'lh.pbt.', 'lh.toroGI*.'));
if isempty(tmp_Toro)
    doToro = false;
else
    doToro        = true;
    numToro       = length(tmp_Toro);
    Toroprefixes  = cell(numToro,1);
    [~, tmp_name] = fileparts(regexprep(strrep(surfaces{1}, 'lh.pbt.', ''), ...
                             {'_T1w_RAS', '_T1w'}, ''));
    for Toro = 1:numToro
        Toroprefixes{Toro,1} = regexprep(tmp_Toro(Toro).name, {tmp_name, ...
                               'lh.toroGI', 'mm.', '_T1w_RAS', '_T1w'}, '');
    end
end

if doToro
    var_names = {'SubjectID';                                              ...
                  'Mean_PBT';                    'SD_PBT';                 ...
                  'Mean_Thickness';              'SD_Thickness';           ...
                  'Mean_Gyrification';           'SD_Gyrification';        ...
                  'Mean_Depth';                  'SD_Depth';               ...
                  'Mean_FractalDimension';       'SD_FractalDimension'};
              
    for Toro = 1:numToro
        var_names{length(var_names)+1,1} = ['Mean_ToroGI', Toroprefixes{Toro}, 'mm'];
        var_names{length(var_names)+1,1} = ['SD_ToroGI',   Toroprefixes{Toro}, 'mm'];
    end
    
    var_names = [var_names;                                             ...
                {'Mean_LH_PBT';               'SD_LH_PBT';              ...
                'Mean_RH_PBT';                'SD_RH_PBT';              ...
                'Mean_LH_Thickness';          'SD_LH_Thickness';        ...
                'Mean_RH_Thickness';          'SD_RH_Thickness';        ...
                'Mean_LH_Gyrification';       'SD_LH_Gyrification';     ...
                'Mean_RH_Gyrification';       'SD_RH_Gyrification';     ...
                'Mean_LH_Depth';              'SD_LH_Depth';            ...
                'Mean_RH_Depth';              'SD_RH_Depth';            ...
                'Mean_LH_FractalDimension';   'SD_LH_FractalDimension'; ...
                'Mean_RH_FractalDimension';   'SD_RH_FractalDimension'}];
            
    for Toro = 1:numToro
        var_names{length(var_names)+1,1} = ['Mean_LH_ToroGI', Toroprefixes{Toro}, 'mm'];
        var_names{length(var_names)+1,1} = ['SD_LH_ToroGI',   Toroprefixes{Toro}, 'mm'];
        var_names{length(var_names)+1,1} = ['Mean_RH_ToroGI', Toroprefixes{Toro}, 'mm'];
        var_names{length(var_names)+1,1} = ['SD_RH_ToroGI',   Toroprefixes{Toro}, 'mm'];        
    end
    
else
    var_names = {'SubjectID';                                             ...
                 'Mean_PBT';                    'SD_PBT';                 ...
                 'Mean_Thickness';              'SD_Thickness';           ...
                 'Mean_Gyrification';           'SD_Gyrification';        ...
                 'Mean_Depth';                  'SD_Depth';               ...
                 'Mean_FractalDimension';       'SD_FractalDimension';    ...
                 'Mean_LH_PBT';                 'SD_LH_PBT';              ...
                 'Mean_RH_PBT';                 'SD_RH_PBT';              ...
                 'Mean_LH_Thickness';           'SD_LH_Thickness';        ...
                 'Mean_RH_Thickness';           'SD_RH_Thickness';        ...
                 'Mean_LH_Gyrification';        'SD_LH_Gyrification';     ...
                 'Mean_RH_Gyrification';        'SD_RH_Gyrification';     ...
                 'Mean_LH_Depth';               'SD_LH_Depth';            ...
                 'Mean_RH_Depth';               'SD_RH_Depth';            ...
                 'Mean_LH_FractalDimension';    'SD_LH_FractalDimension'; ...
                 'Mean_RH_FractalDimension';    'SD_RH_FractalDimension'};
end

% Initialize
globalSurf = cell(num_surfaces, 31+(numToro*6));
globalSurf = cell2table(globalSurf, 'VariableNames', var_names);

%% Compile measures
% Loop over each subject file and summarize values
for subjs  = 1:num_surfaces
        
    % Get subject ID
    [~, tmp_name]               = fileparts(regexprep(strrep(surfaces{subjs}, 'lh.pbt.', ''), {'_T1w_RAS', '_T1w'}, ''));
    globalSurf.SubjectID{subjs} = tmp_name;
    
    % Get PBT info
    dat_LH                        = cat_io_FreeSurfer('read_surf_data', surfaces{subjs});
    dat_RH                        = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'rh.pbt.'));
    globalSurf.Mean_PBT{subjs}    = mean([dat_LH; dat_RH]);
    globalSurf.SD_PBT{subjs}      = std([dat_LH;  dat_RH]);
    globalSurf.Mean_LH_PBT{subjs} = mean(dat_LH);
    globalSurf.SD_LH_PBT{subjs}   = std(dat_LH);
    globalSurf.Mean_RH_PBT{subjs} = mean(dat_RH);
    globalSurf.SD_RH_PBT{subjs}   = std(dat_RH);
    
    % Get cortical thickness
    try
        dat_LH                              = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'lh.thickness.'));
        dat_RH                              = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'rh.thickness.'));
        globalSurf.Mean_Thickness{subjs}    = mean([dat_LH; dat_RH]);
        globalSurf.SD_Thickness{subjs}      = std([dat_LH;  dat_RH]);
        globalSurf.Mean_LH_Thickness{subjs} = mean(dat_LH);
        globalSurf.SD_LH_Thickness{subjs}   = std(dat_LH);
        globalSurf.Mean_RH_Thickness{subjs} = mean(dat_RH);
        globalSurf.SD_RH_Thickness{subjs}   = std(dat_RH);
    catch
        warning(['Cortical thickness file not found for: ', tmp_name]);
    end
    
    % Get gyrification
    try
        dat_LH                                  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'lh.gyrification.'));
        dat_RH                                  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'rh.gyrification.'));
        globalSurf.Mean_Gyrification{subjs}     = mean([dat_LH; dat_RH]);
        globalSurf.SD_Gyrification{subjs}       = std([dat_LH;  dat_RH]);
        globalSurf.Mean_LH_Gyrification{subjs}  = mean(dat_LH);
        globalSurf.SD_LH_Gyrification{subjs}    = std(dat_LH);
        globalSurf.Mean_RH_Gyrification{subjs}  = mean(dat_RH);
        globalSurf.SD_RH_Gyrification{subjs}    = std(dat_RH);
    catch
        warning(['Gyrification file not found for: ', tmp_name]);
    end
    
    % Get sulcal depth
    try
        dat_LH                          = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'lh.depth.'));
        dat_RH                          = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'rh.depth.'));
        globalSurf.Mean_Depth{subjs}    = mean([dat_LH; dat_RH]);
        globalSurf.SD_Depth{subjs}      = std([dat_LH;  dat_RH]);
        globalSurf.Mean_LH_Depth{subjs} = mean(dat_LH);
        globalSurf.SD_LH_Depth{subjs}   = std(dat_LH);
        globalSurf.Mean_RH_Depth{subjs} = mean(dat_RH);
        globalSurf.SD_RH_Depth{subjs}   = std(dat_RH);
    catch
        warning(['Depth file not found for: ', tmp_name]);
    end
    
    % Get fractal dimension
    try
        dat_LH                                      = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'lh.fractaldimension.'));
        dat_RH                                      = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', 'rh.fractaldimension.'));
        globalSurf.Mean_FractalDimension{subjs}     = mean([dat_LH; dat_RH]);
        globalSurf.SD_FractalDimension{subjs}       = std([dat_LH;  dat_RH]);
        globalSurf.Mean_LH_FractalDimension{subjs}  = mean(dat_LH);
        globalSurf.SD_LH_FractalDimension{subjs}    = std(dat_LH);
        globalSurf.Mean_RH_FractalDimension{subjs}  = mean(dat_RH);
        globalSurf.SD_RH_FractalDimension{subjs}    = std(dat_RH);
    catch
        warning(['Depth file not found for: ', tmp_name]);
    end    
    
    % Get Toro's gyrification indices
    if doToro
        for Toro = 1:numToro
            try
                dat_LH  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', ['lh.toroGI', Toroprefixes{Toro}, 'mm.']));
                dat_RH  = cat_io_FreeSurfer('read_surf_data', strrep(surfaces{subjs}, 'lh.pbt.', ['rh.toroGI', Toroprefixes{Toro}, 'mm.']));        
                globalSurf.(['Mean_ToroGI',     Toroprefixes{Toro}, 'mm']){subjs} = mean([dat_LH; dat_RH]);
                globalSurf.(['SD_ToroGI',       Toroprefixes{Toro}, 'mm']){subjs} = std([dat_LH;  dat_RH]);
                globalSurf.(['Mean_LH_ToroGI',  Toroprefixes{Toro}, 'mm']){subjs} = mean(dat_LH);
                globalSurf.(['SD_LH_ToroGI',    Toroprefixes{Toro}, 'mm']){subjs} = std(dat_LH);
                globalSurf.(['Mean_RH_ToroGI',  Toroprefixes{Toro}, 'mm']){subjs} = mean(dat_RH);
                globalSurf.(['SD_RH_ToroGI',    Toroprefixes{Toro}, 'mm']){subjs} = std(dat_RH);
            catch
                warning(['Toro GI ', Toroprefixes{Toro}, ' file not found for: ', tmp_name]);
            end
        end
    end
end

% Write table, if required
if toWrite
    writetable(globalSurf, outname);
end