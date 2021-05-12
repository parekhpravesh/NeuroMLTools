function prep_proportionImages(mwp1, mask, dir_output)
% Function to create proportion of gray matter, white matter, and 
% cerebrospinal fluid images given a set of CAT mwp[1-3] images
%% Inputs
% mwp1:             full path to a location where mwp[1-3]*.nii files
%                   (created by CAT) are saved OR cell type with each row
%                   being a full path to mwp1 file
% mask:             either a NIfTI image that will be used to mask the GM
%                   file (before writing) OR a single number indicating the
%                   threshold to be applied to the image (before
%                   computation of proportion) [optional]
% dir_output:       fullpath to a location where results will be saved
% 
%% Outputs:
% Proportion of GM, WM, and CSF images are written in dir_output
% 
%% Notes:
% For each voxel, proportion of GM is defined as:
% GM./(GM + WM + CSF) at that voxel
% 
% Similarly, proportion of WM and CSF are defined as:
% WM./(GM + WM + CSF) at that voxel
% CSF./(GM + WM + CSF) at that voxel
% 
% Assumes that mwp* for a given subject are all present in the same folder
% 
% Masking is perhaps a good idea because the proportion of GM images might
% "spill" out into CSF regions
% 
% If a mask is explicitly specified, then proportion images are first
% calculated and then masked out; if a threshold is specified, then
% a mask for each image is first calculated using this threshold and then 
% the proportion image is calculated
% 
% Assumes that the mask image (if provided) has the same dimension and
% header as the rest of the images
% 
% Threshold is applied as a less than operation
% 
%% Defaults:
% mask:             ''
% dir_output:       same as first entry in mwp1 or pwd
% 
%% Author(s):
% Parekh, Pravesh
% May 10, 2021
% MBIAL

%% Check inputs
% Check dir_mri
if ~exist('mwp1', 'var') || isempty(mwp1)
    error('Please provide directory having mwp[1-3] files OR list of mwp1 files to work on');
else
    if iscell(mwp1)
        % Build additional file names
        numImages = length(mwp1);
        mwp2 = strrep(mwp1, 'mwp1', 'mwp2');
        mwp3 = strrep(mwp1, 'mwp1', 'mwp3');
        f    = @(x) exist(x, 'file');
        
        % Make sure all mwp[1-3] images exist
        if sum(logical(cellfun(f, mwp1))) ~= numImages
            error('One or more mwp1 files not found');
        end
        if sum(logical(cellfun(f, mwp2))) ~= numImages
            error('One or more mwp2 files not found');
        end
        if sum(logical(cellfun(f, mwp3))) ~= numImages
            error('One or more mwp3 files not found');
        end
    else
        if ~exist(mwp1, 'dir')
            error(['Unable to find: ', mwp1]);
        else
            % Find all NIfTI images
            list_mwp1 = dir(fullfile(mwp1, 'mwp1*.nii'));
            mwp1      = fullfile(mwp1, {list_mwp1(:).name}');
            numImages = length(mwp1);
            
            % Make sure mwp[2-3] images exist
            f    = @(x) exist(x, 'file');
            mwp2 = strrep(mwp1, 'mwp1', 'mwp2');
            mwp3 = strrep(mwp1, 'mwp1', 'mwp3');
            if sum(logical(cellfun(f, mwp2))) ~= numImages
                error('One or more mwp2 files not found');
            end
            if sum(logical(cellfun(f, mwp3))) ~= numImages
                error('One or more mwp3 files not found');
            end
        end
    end
end

% Check mask
if ~exist('mask', 'var') || isempty(mask)
    doMask = false;
else
    doMask = true;
    if isnumeric(mask)
        doThresh = true;
    else
        doThresh = false;
        if ~exist(mask, 'file')
            error(['Unable to find: ', mask]);
        else
            volMask = spm_vol(mask);
            datMask = spm_read_vols(volMask);
        end
    end
end

% Check dir_output
if ~exist('dir_output', 'var') || isempty(dir_output)
    if isempty(fileparts(mwp1{1}))
        dir_output = pwd;
    else
        dir_output = fileparts(mwp1{1});
    end
else
    if ~exist(dir_output, 'dir')
        mkdir(dir_output);
    end
end

%% Loop over each image and compute
for images = 1:numImages
    
    % Get header
    mat_mwp1 = spm_get_space(mwp1{images});
    mat_mwp2 = spm_get_space(mwp2{images});
    mat_mwp3 = spm_get_space(mwp3{images});
    
    % Make sure headers internally match
    if any(mat_mwp1 - mat_mwp2) | any(mat_mwp1 - mat_mwp3) %#ok<OR2>
        error('Header mismatch between mwp[1-3] files');
    else
        if ~doThresh
            % Make sure header matches mask file
            if any(mat_mwp1 - volMask.mat)
                error('Header mismatch between mwp1 and mask file');
            end
        end
    end
    
    % Read data
    vol_mwp1 = spm_vol(mwp1{images});
    dat_mwp1 = spm_read_vols(vol_mwp1);
    dat_mwp2 = spm_read_vols(spm_vol(mwp2{images}));
    dat_mwp3 = spm_read_vols(spm_vol(mwp3{images}));
    
    % Apply threshold, if necessary
    if doThresh
        dat_mwp1(dat_mwp1 < mask) = 0;
        dat_mwp2(dat_mwp2 < mask) = 0;
        dat_mwp3(dat_mwp3 < mask) = 0;
    end
    
    % Calculate proportion image
    denom   = dat_mwp1 + dat_mwp2 + dat_mwp3;
    propGM  = dat_mwp1./denom;
    propWM  = dat_mwp2./denom;
    propCSF = dat_mwp3./denom;
    
    % Apply mask, if necessary
    if doMask
        propGM  = propGM.*datMask;
        propWM  = propWM.*datMask;
        propCSF = propCSF.*datMask;
    end
    
    % Edit header and write out
    [~, tmpName]   = fileparts(mwp1{images});
    vol_mwp1.fname = fullfile(dir_output, ['PropGM_', tmpName, '.nii']);
    spm_write_vol(vol_mwp1, propGM);
    
    vol_mwp1.fname = fullfile(dir_output, ['PropWM_', tmpName, '.nii']);
    spm_write_vol(vol_mwp1, propWM);
    
    vol_mwp1.fname = fullfile(dir_output, ['PropCSF_', tmpName, '.nii']);
    spm_write_vol(vol_mwp1, propCSF);   
end