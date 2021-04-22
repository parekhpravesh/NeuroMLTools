function get_cat_qa(files, dir_output)
% Function to read CAT report XML or .mat files and extract quality 
% assurance measures, percentages, and grades into a csv file
%% Inputs:
% files:            cell with full paths to .XML or .mat files
%                   (rows are filenames) OR full path to a directory
%                   containing these .XML or .mat files
% dir_output:       full path to where QA file is written out
% 
%% Output:
% A csv file QA_CAT_<DDMMMYYYY> is written either in dir_output or the 
% directory where the first .XML/.mat file is; if neither of these 
% locations exist, then it is saved in the pwd
% 
% The QA_CAT_DDMMMYYYY.csv file contains the following:
%   * Subject ID
%   * Voxel resolution
%   * Resolution:
%       - RMS value
%       - RMS RPS
%       - RMS grade
%   * Noise:
%       - NCR value
%       - NCR RPS
%       - NCR grade
%   * Bias:
%       - ICR value
%       - ICR RPS
%       - ICR grade
%   * Weighted IQR:
%       - Weighted IQR value
%       - Weighted IQR RPS
%       - Weighted IQR grade
%   * Surface measures:
%       - Mean surface Euler number
%       - Mean surface defect number
%       - Mean surface defect area
%       - Surface intensity RMSE
%       - Position RMSE
%       - Surface self intersection
%   * Mean tissue intensities:
%       - GM
%       - WM
%       - CSF
%       - Background
%   * Standard deviation tissue intensities:
%       - GM
%       - WM
%       - CSF
%       - Background
%   * Relative mean tissue intensities (relative to maximum of WM/GM)
%       - GM
%       - WM
%       - CSF
%       - Background
%   * Relative standard deviation tissue intensities
%       - GM
%       - WM
%       - CSF
%       - Background
%   * Contrast between tissue classes:
%       - Contrast
%       - Relative contrast
%   * Versions
%       - MATLAB version
%       - SPM version
%       - CAT version
% 
% If XML files are input, CAT automatically saves the XML file as a mat 
% file with the same filename
% 
% If XML files are input and full paths are not provided, MATLAB will
% generate a warning about empty directory names
% 
%% Notes:
% If XML files are input, the function uses cat_io_xml function to read the 
% xml file into a MATLAB structure
% 
% The rest of the code relies on the grades and functions 
% marks2str, mark2rps, and mark2grad defined in cat_main.m file
% (lines 2051-2054; CAT 12.5 v1363)
% 
% These are also defined in cat_stat_marks and the XML file mentions
% cat_vol_qa functions for quality assurance
% 
% Resolution RMS:   RMS error of voxel size
% Noise NCR:        Noise to contrast ratio
% Bias ICR:         Inhomogeneity to contrast ratio
% 
% value:            The value passed to mark2rps
% rps:              Percentage rating points
% grade:            A+ to F
% 
% See S.ratings_help.qualitymeasures for a description of some of these 
% quality measures
% 
% Additionally, see cat_vol_qa for calculation details, especially
% pertaining to tissue intensities and contrast
% 
%% Default:
% dir_output:      same as first input file OR pwd
% 
%% Author(s):
% Parekh, Pravesh
% January 03, 2019
% April 22, 2021
% MBIAL

%% Check inputs and assign defaults
% Check files
if ~exist('files', 'var') || isempty(files)
    error('files should be provided');
else
    if iscell(files)
        num_files = length(files);
        % Make sure all files exist
        f = @(x) exist(x, 'file');
        if sum(logical(cellfun(f, files))) ~= num_files
            error('One or more files not found');
        end
    else
        if ~exist(files, 'dir')
            error(['Unable to find directory: ', files]);
        else
            tmp_files = dir(fullfile(files, 'cat_*.mat'));
            tmp_files = {tmp_files(:).name}';
            files     = fullfile(files, tmp_files);
            num_files = length(files);
        end
    end
end

% Check dir_output
if ~exist('dir_output', 'var') || isempty(dir_output)
    % Check if full path is present in the first file
    if isempty(fileparts(files{1}))
        outname = fullfile(pwd, ['QA_CAT_', datestr(now, 'ddmmmyyyy'), '.csv']);
    else
        outname = fullfile(fileparts(files{1}), ['QA_CAT_', datestr(now, 'ddmmmyyyy'), '.csv']); 
    end
else
    if ~exist(dir_output, 'dir')
        mkdir(dir_output);
    end
    outname = fullfile(dir_output, ['QA_CAT_', datestr(now, 'ddmmmyyyy'), '.csv']);
end

%% Define functions
% These lines are taken from cat_main.m (lines 2051-2054; CAT 12.5 v1363)
% These lines can now be found in cat_main_reportcmd.m (lines 23-25;  
% 1701 2020-08-24) and in other places
grades      = {'A+','A','A-','B+','B','B-','C+','C','C-','D+','D','D-','E+','E','E-','F'};
mark2rps    = @(mark) min(100,max(0,105 - mark*10)) + isnan(mark).*mark;
mark2grad   = @(mark) grades{min(numel(grades),max(max(isnan(mark)*numel(grades),1),round((mark+2/3)*3-3)))};
  
%% Initialize
header = {'SubjectID';           'voxelResolution';                               ...
          'resolutionRMS_value'; 'resolutionRMS_rps';    'resolutionRMS_grade';   ...
          'noiseNCR_value';      'noiseNCR_rps';         'noiseNCR_grade';        ...
          'biasICR_value';       'biasICR_rps';          'biasICR_grade';         ...
          'weightedIQR_value';   'weightedIQR_rps';      'weightedIQR_grade';     ...
          'surf_EulerNumber';    'surf_defectNumber';    'surf_defectArea';       ...
          'surf_intensityRMSE';  'surf_positionRMSE';    'surf_selfIntersection'; ...
          'meanIntensity_GM';    'meanIntensity_WM';     'meanIntensity_CSF';     ...
          'meanIntensity_BG';    'SDIntensity_GM';       'SDIntensity_WM';        ...
          'SDIntensity_CSF';     'SDIntensity_BG';       'relMeanIntensity_GM';   ...
          'relMeanIntensity_WM'; 'relMeanIntensity_CSF'; 'relMeanIntensity_BG';   ...
          'relSDIntensity_GM';   'relSDIntensity_WM';    'relSDIntensity_CSF';    ...
          'relSDIntensity_BG';   'Contrast';             'relativeContrast'; 
          'versionMATLAB';       'versionSPM';           'versionCAT'};
      
measures = cell(num_files, length(header));

%% Get values
for file = 1:num_files
    
    % Check if XML or mat file
    ext  = strsplit(files{file}, '.');
    ext  = ext{end};
    if strcmpi(ext, 'xml')
        str  = cat_io_xml(files{file});
    else
        str = load(files{file}, 'S');
        str = str.S;
    end
        
    % Subject ID
    [~, measures{file,1}]  = fileparts(regexprep(files{file}, {'cat_', '_T1w_RAS', '_T1w', '.mat'}, ''));
    
    % Skip subject if there was an error
    if isfield(str, 'error')
        warning(['Error in file: ', files{file}, '; skipping this file']);
        continue;
    end    
    
    % Voxel resolution
    measures{file,2}  = [num2str(str.qualitymeasures.res_vx_vol(1), '%.4f'), ' × ', ...
                         num2str(str.qualitymeasures.res_vx_vol(2), '%.4f'), ' × ', ...
                         num2str(str.qualitymeasures.res_vx_vol(3), '%.4f')];    
    
    % Resolution
    measures{file,3}  = str.qualityratings.res_RMS;
    measures{file,4}  = mark2rps(str.qualityratings.res_RMS);
    measures{file,5}  = mark2grad(str.qualityratings.res_RMS);
    
    % Noise
    measures{file,6}  = str.qualityratings.NCR;
    measures{file,7}  = mark2rps(str.qualityratings.NCR);
    measures{file,8}  = mark2grad(str.qualityratings.NCR);
    
    % Bias
    measures{file,9}  = str.qualityratings.ICR;
    measures{file,10} = mark2rps(str.qualityratings.ICR);
    measures{file,11} = mark2grad(str.qualityratings.ICR);
    
    % Weighted IQR
    measures{file,12} = str.qualityratings.IQR;
    measures{file,13} = mark2rps(str.qualityratings.IQR);
    measures{file,14} = mark2grad(str.qualityratings.IQR);
    
    try
        % Surface Euler number
        measures{file,15} = str.qualitymeasures.SurfaceEulerNumber;
        
        % Surface defect number
        measures{file,16} = str.qualitymeasures.SurfaceDefectNumber;
        
        % Surface defect area
        measures{file,17} = str.qualitymeasures.SurfaceDefectArea;
        
        % Surface intensity RMSE
        measures{file,18} = str.qualitymeasures.SurfaceIntensityRMSE;
        
        % Surface position RMSE
        measures{file,19} = str.qualitymeasures.SurfacePositionRMSE;
        
        % Surface self intersection
        measures{file,20} = str.qualitymeasures.SurfaceSelfIntersections;
        
    catch
        % Set surface measures to NaN
        measures{file,15} = NaN;        
        measures{file,16} = NaN;
        measures{file,17} = NaN;
        measures{file,18} = NaN;
        measures{file,19} = NaN;
        measures{file,20} = NaN;
    end
    
    % Mean tissue intensities
    measures{file,21} = str.qualitymeasures.tissue_mn(3);
    measures{file,22} = str.qualitymeasures.tissue_mn(4);
    measures{file,23} = str.qualitymeasures.tissue_mn(2);
    measures{file,24} = str.qualitymeasures.tissue_mn(1);
    
    % Standard deviation tissue intensities
    measures{file,25} = str.qualitymeasures.tissue_std(3);
    measures{file,26} = str.qualitymeasures.tissue_std(4);
    measures{file,27} = str.qualitymeasures.tissue_std(2);
    measures{file,28} = str.qualitymeasures.tissue_std(1);
    
    % Mean relative tissue intensities
    measures{file,29} = str.qualitymeasures.tissue_mnr(3);
    measures{file,30} = str.qualitymeasures.tissue_mnr(4);
    measures{file,31} = str.qualitymeasures.tissue_mnr(2);
    measures{file,32} = str.qualitymeasures.tissue_mnr(1);
    
    % Standard deviation relative tissue intensities
    measures{file,33} = str.qualitymeasures.tissue_stdr(3);
    measures{file,34} = str.qualitymeasures.tissue_stdr(4);
    measures{file,35} = str.qualitymeasures.tissue_stdr(2);
    measures{file,36} = str.qualitymeasures.tissue_stdr(1);
    
    % Contrast between tissue classes
    measures{file,37} = str.qualitymeasures.contrast;
    measures{file,38} = str.qualitymeasures.contrastr;
    
    % MATLAB version
    measures{file,39} = num2str(str.software.version_matlab);
    
    % SPM version
    measures{file,40} = num2str(str.software.version_spm);
    
    % CAT version
    measures{file,41} = [num2str(str.software.version_cat), ' v', ...
                         num2str(str.software.revision_cat)];
                     
    % Clear up for the next file
    clear str
end

%% Convert to table and save as csv file
measures = cell2table(measures, 'VariableNames', header);
writetable(measures, outname);