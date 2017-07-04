warning('off','all');
%% project directory
PROJpath = '/Volumes/JetDrive/data/SJ_Ahn/BM/';

% target directory to save NIFTI file
DATApath = fullfile(PROJpath,'metastasis');


% specify list
fn_list = fullfile(PROJpath,'Breast cancer mets 20170121.xlsx');
T = readtable(fn_list);

subjlist = T.serialNumber;

% Independent variable for LR
%--------------------------------------------------------------------------
Age = T.age;
luminal = T.luminal;
HER2 = T.HER2;
TripleNega=T.Basal;
X = [Age, TripleNega,HER2,luminal];
geneType={'TripleNega','HER2','luminal'};


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
PROB = zeros(nvox,nsubj);
BMidx = zeros(nvox,1);
parfor i=1:nvox,
    Y = IMG(:,idbrainmask(i))>0;
    if sum(Y)<1, continue; end
    prob = logisticRegress(X,Y);
    PROB(i,:) = prob;
    BMidx(i) = 1;
end



% Set output file name (Statistical significance)
%--------------------------------------------------------------------------
OUTpath = fullfile(PROJpath,'analysis','LR'); mkdir(OUTpath);

for i=1:length(geneType),
    filename = sprintf('logistic_log10_P_uncorr_%s.nii',geneType{i});
    fn_out  = fullfile(OUTpath,filename);
    
    % Write resulting overlap imaging (number of overlaps)
    Iout       = zeros(vref.dim);
    Iout(idbrainmask) = -log10(Pval);
    vout       = vo;
    vout.fname = fn_out;
    vout.dt    = [16 0];
    spm_write_vol(vout,Iout);
end
