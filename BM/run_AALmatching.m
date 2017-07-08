warning('off','all');
% project directory
PROJpath = '/Volumes/JetDrive/data/SJ_Ahn/BM/';

% target directory to save NIFTI file
DATApath = fullfile(PROJpath,'metastasis');

% Set radius
RADIUS = 15;

% specify list
%--------------------------------------------------------------------------
fn_list = fullfile(PROJpath,'Breast cancer mets 20170121.xlsx');
T = readtable(fn_list);
subjlist = T.serialNumber;
Age = T.age;
Dx  = T.Dx;


% Load template
%--------------------------------------------------------------------------
fn_atlas = '/Volumes/JetDrive/data/SJ_Ahn/matlabscripts/atlas/AAL_157x189x156.nii';
% fn_atlas = '/Volumes/JetDrive/data/SJ_Ahn/matlabscripts/atlas/BM_atlas_final_157x189x156.nii';
% fn_atlas = '/Volumes/JetDrive/data/SJ_Ahn/matlabscripts/atlas/hammers_157x189x156.nii';
vo_atlas = spm_vol(fn_atlas);
atlas = spm_read_vols(vo_atlas); nroi = max(atlas(:));


%
nsubj = length(subjlist);
fn_ref = fullfile(DATApath,num2str(subjlist(1)),['wmetastasis_roi_' num2str(RADIUS) 'mm.nii']);
vref = spm_vol(fn_ref);
[idbrainmask,idgm,idwm,idcsf] = fmri_load_maskindex(vref);

idremove = [];
CNT = zeros(nsubj,nroi);
cnt = 1;
for c=1:nsubj,
    
    % Name of participants
    %----------------------------------------------------------------------
    subjname = num2str(subjlist(c));
    fprintf('[%03d/%03d], subj-%s is running...\n',c,nsubj,subjname);
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
    
    for i=1:nroi,
        idx = find(atlas==i);
        CNT(c,i) = sum(ROI(idx));
    end
end
Dx(idremove) = [];
CNT(idremove,:) = [];



% Triple Negative: Crosstab analysis
%--------------------------------------------------------------------------
vx_thr = 0.2*pi*RADIUS*RADIUS;
Pval_TN = zeros(nroi,1);
Chi2_TN = zeros(nroi,1);
fn = 'basal.csv';
fid = fopen(fn,'w+');
fprintf(fid,'basal (pct),non-basal (pct), chi2, p\n');
for i=1:nroi,
    BM_at_roi = CNT(:,i)>vx_thr;
    [tbl,chi2,p] = crosstab(BM_at_roi,T.Basal);
    Pval_TN(i) = p;
    Chi2_TN(i) = chi2;
    
    fprintf(fid,'%.2f, %.2f, %.2f, %.3f\n',100*mean(BM_at_roi(T.Basal==1)) ,100*mean(BM_at_roi(T.Basal==0)),chi2,p);
end
fclose(fid);
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(Pval_TN);
[find(Pval_TN<0.05) Pval_TN(Pval_TN<0.05) adj_p(Pval_TN<0.05)];
length(find(Pval_TN<0.05))
PCT_Basal = 100*sum(CNT(T.Basal==1,:)>vx_thr)/sum(T.Basal);
[val,id]=sort(PCT_Basal,'descend');



% HER2+: Crosstab analysis
%--------------------------------------------------------------------------
Pval_HER2 = zeros(nroi,1);
Chi2_HER2 = zeros(nroi,1);
fn = 'her2.csv';
fid = fopen(fn,'w+');
fprintf(fid,'her2 (pct),non-her2 (pct), chi2, p\n');
for i=1:nroi,
    BM_at_roi = CNT(:,i)>vx_thr;
    [tbl,chi2,p] = crosstab(BM_at_roi,T.HER2);
    Pval_HER2(i) = p;
    Chi2_HER2(i) = chi2;
    
    fprintf(fid,'%.2f, %.2f, %.2f, %.3f\n',100*mean(BM_at_roi(T.HER2==1)) ,100*mean(BM_at_roi(T.HER2==0)),chi2,p);
end
fclose(fid);
% [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(Pval_HER2);
adj_p = Pval_HER2*116;
[find(Pval_HER2<0.05) Pval_HER2(Pval_HER2<0.05) adj_p(Pval_HER2<0.05)];
length(find(Pval_HER2<0.05))
PCT_HER2 = 100*sum(CNT(T.HER2==1,:)>vx_thr)/sum(T.HER2);
[val,id]=sort(PCT_HER2,'descend');


% Luminal: Crosstab analysis
%--------------------------------------------------------------------------
Pval_LU = zeros(nroi,1);
Chi2_LU = zeros(nroi,1);
fn = 'luminal.csv';
fid = fopen(fn,'w+');
fprintf(fid,'luminal (pct),non-luminal (pct), chi2, p\n');
for i=1:nroi,
    BM_at_roi = CNT(:,i)>vx_thr;
    [tbl,chi2,p] = crosstab(BM_at_roi,T.luminal);
    Pval_LU(i) = p;
    Chi2_LU(i) = chi2;
    
    fprintf(fid,'%.2f, %.2f, %.2f, %.3f\n',100*mean(BM_at_roi(T.luminal==1)) ,100*mean(BM_at_roi(T.luminal==0)),chi2,p);
end
fclose(fid);
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(Pval_LU);
[find(Pval_LU<0.05) Pval_LU(Pval_LU<0.05) adj_p(Pval_LU<0.05)];
length(find(Pval_LU<0.05))
PCT_Luminal = 100*sum(CNT(T.luminal==1,:)>vx_thr)/sum(T.luminal);
[val,id]=sort(PCT_Luminal,'descend');




return

vo_out = vo_atlas;
vo_out.fname='overlap_all.nii';
vo_out.dt=[16,0];
IMG = zeros(vo_atlas.dim);
pctOverlap = 100*mean(CNT>vx_thr);
for i=1:nroi,
    idx = find(atlas==i);
    IMG(idx)=pctOverlap(i);
end
spm_write_vol(vo_out,IMG);

vo_out = vo_atlas;
vo_out.fname='overlap_basal.nii';
vo_out.dt=[16,0];
IMG = zeros(vo_atlas.dim);
pctOverlap = 100*mean(CNT(T.Basal==1,:)>vx_thr);
for i=1:nroi,
    idx = find(atlas==i);
    IMG(idx)=pctOverlap(i);
end
spm_write_vol(vo_out,IMG);


vo_out = vo_atlas;
vo_out.fname='overlap_her2.nii';
vo_out.dt=[16,0];
IMG = zeros(vo_atlas.dim);
pctOverlap = 100*mean(CNT(T.HER2==1,:)>vx_thr);
for i=1:nroi,
    idx = find(atlas==i);
    IMG(idx)=pctOverlap(i);
end
spm_write_vol(vo_out,IMG);

vo_out = vo_atlas;
vo_out.fname='overlap_luminal.nii';
vo_out.dt=[16,0];
IMG = zeros(vo_atlas.dim);
pctOverlap = 100*mean(CNT(T.luminal==1,:)>vx_thr);
for i=1:nroi,
    idx = find(atlas==i);
    IMG(idx)=pctOverlap(i);
end
spm_write_vol(vo_out,IMG);