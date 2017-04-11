
%% project directory
% PROJpath = '/home/ahn/metastasis/';
PROJpath = '/Volumes/JetDrive/data/SJ_Ahn/';

% target directory to save NIFTI file
DATApath = fullfile(PROJpath,'metastasis');

% Set radius of BM
RADIUS = 15;

% specify list
fn_list = fullfile(PROJpath,'Breast cancer mets 20170121.xlsx');
% fn_list = fullfile(PROJpath,'BM_basal-.xlsx');
[a,b,xlsData] = xlsread(fn_list);
subjlist = cell2mat(xlsData(2:end,3));


%%
%--------------------------------------------------------------------------
% DO NOT CHANGE BELOW
%--------------------------------------------------------------------------
%%
nsubj = length(subjlist);
fname = sprintf('wmetastasis_roi_%dmm.nii',RADIUS);
fn_ref = fullfile(DATApath,num2str(subjlist(1)),fname);
vref = spm_vol(fn_ref);
[idbrainmask,idgm,idwm,idcsf] = fmri_load_maskindex(vref);

Overlap = zeros(vref.dim);
cnt = 0;
for c=1:nsubj,
    
    % Name of participants
    %----------------------------------------------------------------------
    subjname = num2str(subjlist(c));
    
    % Load ROI image
    %----------------------------------------------------------------------
    fn = fullfile(DATApath,subjname,fname);
    if ~exist(fn,'file'),
        fprintf('%s, is not included.\n',subjname);
        continue;
    end;
    vo  = spm_vol(fn);
    ROI = spm_read_vols(vo);
    ROI(ROI>0) = 1;
    
    % Load ROI image
    %----------------------------------------------------------------------
    Overlap(idbrainmask) = Overlap(idbrainmask) + ROI(idbrainmask);
    cnt = cnt+1;
end


% Writing results
%--------------------------------------------------------------------------

OUTpath = fullfile(PROJpath,'analysis'); mkdir(OUTpath);

% Write resulting overlap imaging (number of overlaps)
fn_out     = fullfile(OUTpath,'overlap_count.nii');
vout       = vo;
vout.fname = fn_out;
vout.dt    = [2 0];
spm_write_vol(vout,Overlap);


% Write resulting overlap imaging (number of overlaps)
fn_out     = fullfile(OUTpath,'overlap_percent.nii');
IMG        = zeros(vo.dim);
vout       = vo;
vout.fname = fn_out;
vout.dt    = [2 0];
spm_write_vol(vout,100*Overlap/cnt);
