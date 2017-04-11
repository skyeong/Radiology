function [idbrainmask,idgm,idwm,idcsf] = fmri_load_maskindex(vref)

% journal: DPARSF: A MATLAB Toolbox for Pipeline Data Analysis of Resting-State fMRI
% nearest neighbour
idbrainmask = BM_resample_from(vref,'mask_ICV.nii');  
idgm        = BM_resample_from(vref,'GreyMask_02_91x109x91.img');
idwm        = BM_resample_from(vref,'WhiteMask_09_91x109x91.img');
idcsf       = BM_resample_from(vref,'CsfMask_07_91x109x91.img');


function idroi = BM_resample_from(vref,mask_name)
BMpath = fileparts(which('fmri_load_maskindex.m'));

if nargin<2,
    error('Modality should be correctly specified ....');
end

% Read reference image
if ischar(vref),Masking
    vo_ref = spm_vol(vref);
else
    vo_ref = vref;
end;

% Mask image: Get XYZ-Coordinates
fn_mask = fullfile(BMpath,'brainmask',mask_name);
vo_mask = spm_vol(fn_mask);
MASK    = spm_read_vols(vo_mask);
idmask  = find(MASK>0);
[vx, vy, vz] = ind2sub(vo_mask.dim,idmask);
Vxyz = [vx, vy, vz, ones(size(vx,1),1)];
Rxyz = (vo_mask.mat)*Vxyz';

% Resampled in standard normalized space
Vxyz = round(pinv(vo_ref.mat)*Rxyz);

% Transformation from original MASK space to the reference space
idroi = sub2ind(vo_ref.dim, Vxyz(1,:), Vxyz(2,:), Vxyz(3,:));
idroi = unique(idroi);
