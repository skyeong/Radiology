function dcm2nii(workdir, fn_dcm2nii)
% dcm2nii (workdir, dcm2nii)

% Change working directory
prevdir = pwd;
cd(workdir);

% converting images and remove unnecessary images
cmd = sprintf('! %s -g n -4 *.dcm',fn_dcm2nii); eval(cmd);
eval('!rm -rf co*.nii o*.nii');
eval('!mv *.nii metastasis.nii');

% Return to previous directory
cd(prevdir)
