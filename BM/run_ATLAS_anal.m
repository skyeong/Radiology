
%% project directory
% PROJpath = '/home/ahn/metastasis/';
PROJpath = '/Volumes/JetDrive/data/SJ_Ahn/';

% target directory to save NIFTI file
DATApath = fullfile(PROJpath,'metastasis');

% specify list
fn_list = fullfile(PROJpath,'Breast cancer mets 20170106.xlsx');
[a,b,xlsData] = xlsread(fn_list);
subjlist = cell2mat(xlsData(2:end,3));


%%
%--------------------------------------------------------------------------
% DO NOT CHANGE BELOW
%--------------------------------------------------------------------------
%%
nsubj = length(subjlist);
fn_ref = fullfile(DATApath,num2str(subjlist(1)),'wmetastasis_roi_10mm.nii');

% ATLAS
fn_atlas = fullfile(PROJpath,'atlas','rBM_atlas_final.nii');
% fn_atlas = fullfile(PROJpath,'atlas','rBM_atlas_final.nii');
vo_atlas = spm_vol(fn_atlas);
ATLAS = spm_read_vols(vo_atlas);
nrois = max(ATLAS(:));

% Check existance of BM for each ROI
BM = zeros(nsubj,nrois);
for c=1:nsubj,
     % Name of participants
    %----------------------------------------------------------------------
    subjname = num2str(subjlist(c));
    
    % Load ROI image
    %----------------------------------------------------------------------
    fn = fullfile(DATApath,subjname,'wmetastasis_roi_10mm.nii');
    if ~exist(fn,'file'),
        fprintf('%s, is not included.\n',subjname);
        continue;
    end;
    vo  = spm_vol(fn);
    IMG = spm_read_vols(vo);

    % Split cluster into individual ROIs
    idloc = find(IMG>0);
    [vx, vy, vz] = ind2sub(vo.dim, idloc);
    ilab = spm_clusters([vx, vy, vz]')';

    % Iterate for each ROI
    for i=1:max(ilab),
        idx = find(ilab==i);
        idroi = idloc(idx);

        overlap = zeros(nrois,1);
        for j=1:nrois,
            idarea = find(ATLAS==j);
            idcheck = intersect(idarea,idroi);
            overlap(j) = length(idcheck)/length(idroi);
        end
        [val,id] = max(overlap);
        BM(c,id) = BM(c,id) + 1;
    end
end
