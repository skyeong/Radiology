
%% project directory
% PROJpath = '/home/ahn/metastasis';
PROJpath = '/Volumes/JetDrive/data/SJ_Ahn/';

% target directory to save NIFTI file
DATApath = fullfile(PROJpath,'metastasis');

% location of excel file
LOCpath  = fullfile(PROJpath,'coordinates');

% dcm2nii fullpath
dcm2nii_exe = '/Applications/mricron/dcm2nii';

% Set radius
RADIUS = 15;

% specify list
fn_list = fullfile(PROJpath,'Breast cancer mets 20170106.xlsx');
[a,b,xlsData] = xlsread(fn_list);
subjlist = cell2mat(xlsData(2:end,3));





%%
%--------------------------------------------------------------------------
% DO NOT CHANGE BELOW
%--------------------------------------------------------------------------

%%
warning('off','all');
nsubj = length(subjlist);
for c=1:nsubj,
    
    % Name of participants
    %----------------------------------------------------------------------
    subjname = num2str(subjlist(c));
    fprintf('[%03d/%03d] %s is loaded ...\n',c,nsubj,subjname);

    % Convert dicom to nifti
    %----------------------------------------------------------------------
    path_source = fullfile(PROJpath,'dicom',subjname);
    %dcm2nii(path_source, dcm2nii_exe);
    
    % Change name of file
    %----------------------------------------------------------------------
    path_target = fullfile(PROJpath,'metastasis',subjname); mkdir(path_target);
    fn_source = fullfile(path_source,'metastasis.nii');
    fn_target = fullfile(path_target,'metastasis.nii');
    %status = copyfile(fn_source,fn_target,'f');
    %if status,
    %   fprintf('  : dicom image was converted.\n');
    %   delete(fn_source);
    %end
  
    
    % Load image
    %----------------------------------------------------------------------
    fn = fullfile(DATApath,subjname,'metastasis.nii');
    vo = spm_vol(fn);
    DIM = vo.dim;
    ROI = zeros(DIM);
    
    
    % Load center coordinates of metastasis
    %----------------------------------------------------------------------
    fprintf('  : create individual ROI maps using center coordinates.\n');
    fn_xls = fullfile(LOCpath,sprintf('%s.xls',subjname));
    [a,b,xlsData] = xlsread(fn_xls);
    Vxyz = round(cell2mat(xlsData(2:end,3:5)));
    nroi = size(Vxyz,1);
    for i=1:nroi,
        x = Vxyz(i,1);
        y = DIM(2)-Vxyz(i,2);
        z = Vxyz(i,3);
        ROI(x-2:x+2,y-2:y+2,z-2:z+2)=1;
    end
    
    
    % Write roi center coordinate into image
    %----------------------------------------------------------------------
    fn_roi = fullfile(DATApath,subjname,'metastasis_roi.nii');
    vout = vo;
    vout.fname = fn_roi;
    spm_write_vol(vout,ROI);
    
    
    % Normalize roi image
    %----------------------------------------------------------------------
    fprintf('  : normalizing metastasis.nii and metastasis_roi.nii.\n');
    normalization(fn_target, fn_roi)
    
    
    % Draw 10-mm radius ROI in MNI space
    %----------------------------------------------------------------------
    fprintf('  : drawing %d-mm sphere ROI.\n',RADIUS);
    fn_wroi = fullfile(DATApath,subjname,'wmetastasis_roi.nii');
    vo_wroi = spm_vol(fn_wroi);
    ROI = spm_read_vols(vo_wroi);
    idroi = find(ROI>0);
    nroi2 = length(idroi);
    
    IMG = zeros(vo_wroi.dim);
    for i=1:nroi2
        [vx, vy, vz] = ind2sub(vo_wroi.dim, idroi(i));
        Vxyz = [vx vy vz 1];
        Rxyz = vo_wroi.mat*Vxyz';
        idmask = generate_rois(Rxyz(1:3)', RADIUS, fn_wroi);
        IMG(idmask) = 1;
    end
    
    
    % Write image
    %----------------------------------------------------------------------
    fname = sprintf('wmetastasis_roi_%dmm.nii',RADIUS);
    fn_out = fullfile(DATApath,subjname,fname);
    vout = vo_wroi;
    vout.fname=fn_out;
    spm_write_vol(vout,IMG);
    
end