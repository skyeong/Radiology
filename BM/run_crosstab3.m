
%% project directory
PROJpath = '/Volumes/JetDrive/data/SJ_Ahn/BM/';

% target directory to save NIFTI file
DATApath = fullfile(PROJpath,'metastasis');

% Set radius
RADIUS = 15;

% specify list
fn_list = fullfile(PROJpath,'Breast cancer mets 20170121.xlsx');
[a,b,xlsData] = xlsread(fn_list);
hdrs = xlsData(1,:);
dats = xlsData(2:end,:);

subjlist = cell2mat(dats(:,find_column_number(hdrs,'serial number')));
Dx = cell2mat(dats(:,find_column_number(hdrs,'Dx')));



%%
%--------------------------------------------------------------------------
% DO NOT CHANGE BELOW
%--------------------------------------------------------------------------

%%
nsubj = length(subjlist);
fn_ref = fullfile(DATApath,num2str(subjlist(1)),'wmetastasis_roi_10mm.nii');
vref = spm_vol(fn_ref);
[idbrainmask,idgm,idwm,idcsf] = fmri_load_maskindex(vref);

idremove = [];
IMG = zeros([nsubj,prod(vref.dim)]);
cnt = 1;
for c=1:nsubj,
    
    % Name of participants
    %----------------------------------------------------------------------
    subjname = num2str(subjlist(c));
    
    
    % Load ROI image
    %----------------------------------------------------------------------
    filename = sprintf('wmetastasis_roi_%dmm.nii',RADIUS);
    fn = fullfile(DATApath,subjname,filename);
    if ~exist(fn,'file'),
        fprintf('%s, is not included.\n',subjname);
        idremove = [idremove; c];
        continue;
    end;
    vo = spm_vol(fn);
    ROI = spm_read_vols(vo);
    
    
    % Load ROI image
    %----------------------------------------------------------------------
    IMG(c,idbrainmask) = ROI(idbrainmask);
end
Dx(idremove) = [];
IMG(idremove,:) = [];


% Crosstab analysis
%--------------------------------------------------------------------------
nvox = length(idbrainmask);
Pval = zeros(nvox,1);
Chi2 = zeros(nvox,1);
parfor i=1:nvox,
    idx = idbrainmask(i);
    BM_at_voxel  = IMG(:,idx);
    [tbl,chi2,p] = crosstab(BM_at_voxel,Dx);
    Pval(i) = p;
    Chi2(i) = chi2;
end

% Set output file name (number of overlaps)
%--------------------------------------------------------------------------
OUTpath = fullfile(PROJpath,'analysis','3grp'); mkdir(OUTpath);
fn_out  = fullfile(OUTpath,'crosstab_log10_P_uncorr.nii');

% Write resulting overlap imaging (number of overlaps)
Iout       = zeros(vref.dim);
Iout(idbrainmask) = -log10(Pval);
vout       = vo;
vout.fname = fn_out;
vout.dt    = [16 0];
spm_write_vol(vout,Iout);


% Set output file name (percent scale)
%--------------------------------------------------------------------------
OUTpath = fullfile(PROJpath,'analysis','3grp'); mkdir(OUTpath);
fn_out  = fullfile(OUTpath,'crosstab_chi2.nii');


% Write resulting overlap imaging (number of overlaps)
Iout       = zeros(vref.dim);
Iout(idbrainmask) = Chi2;
vout       = vo;
vout.fname = fn_out;
vout.dt    = [16 0];
spm_write_vol(vout,Iout);
