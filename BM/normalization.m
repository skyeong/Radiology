function normalization(fn_metastasis, fn_roi)

fns_resample = '';
fns_resample{1} = fn_metastasis;
fns_resample{2} = fn_roi;

clear matlabbatch;
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {fn_metastasis};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = [fns_resample'];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {fullfile(spm('dir'),'tpm','TPM.nii')};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'eastern';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                             78 76 85];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 0;  % nearest neighbor

% spm_jobman('interactive',matlabbatch);
spm_jobman('run',matlabbatch);