
%% project directory
% PROJpath = '/home/ahn/metastasis/';
PROJpath = '/Volumes/JetDrive/data/SJ_Ahn/';

% target directory to save NIFTI file
DATApath = fullfile(PROJpath,'metastasis');

% ROI size
RADIUS = 15;


% specify list
fn_list = fullfile(PROJpath,'Breast cancer mets 20170121.xlsx');
[a,b,xlsData] = xlsread(fn_list);
hdrs = xlsData(1,:);
dats = xlsData(2:end,:);

subjlist = cell2mat(dats(:,find_column_number(hdrs,'serial number')));
HER2 = cell2mat(dats(:,find_column_number(hdrs,'HER2')));
Luminal = cell2mat(dats(:,find_column_number(hdrs,'luminal')));
Basal = cell2mat(dats(:,find_column_number(hdrs,'Basal')));
CL = Luminal;

% ROIcenters = [-31	-71	-38; 4	-42	-23; 44	18	12; -34	-79	13; -27	2	38;-16	-59	63];


%%
%--------------------------------------------------------------------------
% DO NOT CHANGE BELOW
%--------------------------------------------------------------------------
%%
nsubj = length(subjlist);
fname = sprintf('wmetastasis_roi_%dmm.nii',RADIUS);
fn_ref = fullfile(DATApath,num2str(subjlist(1)),fname);

% Make ROIs
% ROIdir = fullfile(PROJpath,'ROIs','Hippocampus'); mkdir(ROIdir);
ROIdir = fullfile(PROJpath,'ROIs','Luminal'); mkdir(ROIdir);
nrois = 1;
% nrois = length(ROIcenters);
% for i=1:nrois,
%     xyz = ROIcenters(i,:);
%     fn_roi = sprintf('ROI_%d_%d_%d.nii',xyz);
%     fn_out = fullfile(ROIdir,fn_roi);
%     create_ROI_mask(xyz, 3, fn_out, fn_ref,'sphere');
% end

% Check existance of BM for each ROI
BM = zeros(nsubj,nrois);
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
    IMG = spm_read_vols(vo);
    idbm = find(IMG>0);
    
    for i=1:nrois,
        %xyz = ROIcenters(i,:);
        %fn_roi = sprintf('rAAL_%03d.nii',i+36);
        fn_roi = sprintf('ROI%02d.nii',i);
        fn_roi = fullfile(ROIdir,fn_roi);
        vo_roi = spm_vol(fn_roi);
        ROI = spm_read_vols(vo_roi);
        idroi = find(ROI>0);
        idcheck = intersect(idbm,idroi);
        BM(c,i) = length(idcheck)/length(idroi);
    end
end
BM(BM>0)=1;

dat0 = BM(CL==0,:);
dat1 = BM(CL==1,:);

fprintf('CL0, CL1\n');
for i=1:nrois,
   fprintf('%.1f, %.1f\n',100*mean(dat0(:,i)), 100*mean(dat1(:,i))) ;
end
% dlmwrite('~/Desktop/BM.csv',BM)