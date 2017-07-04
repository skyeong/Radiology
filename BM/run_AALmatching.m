warning('off','all');
% project directory
PROJpath = '/Volumes/JetDrive/data/SJ_Ahn/BM/';

% target directory to save NIFTI file
DATApath = fullfile(PROJpath,'metastasis');

% Set radius
RADIUS = 15;

% specify list
fn_list = fullfile(PROJpath,'Breast cancer mets 20170121.xlsx');
T = readtable(fn_list);

subjlist = T.serialNumber;

% Independent variable for LR
%--------------------------------------------------------------------------
Age = T.age;
Dx  = T.Dx;
% X = [Age, TripleNega,HER2,luminal];
geneType={'TripleNega','HER2','luminal'};
Y = cell(0);
for i=1:100,
    if  T.Basal(i)==1,
        Y{i} = 'TripleNega';
    elseif T.HER2(i) ==1,
        Y{i} = 'HER2';
    elseif T.luminal(i)==1,
        Y{i} = 'luminal';
    end
end
Y = categorical(Y)';


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


% Crosstab analysis (3groups)
%--------------------------------------------------------------------------
% Dx: 1(luminal), 2(basal), 3(her2)
unitVol = 0.1*RADIUS*RADIUS*pi;
Pval = zeros(nroi,1);
Chi2 = zeros(nroi,1);
for i=1:nroi,
    BM_at_roi = CNT(:,i)>unitVol;
    [tbl,chi2,p] = crosstab(BM_at_roi,Dx);
    Pval(i) = p;
    Chi2(i) = chi2;
    
    if p<0.05,
       fprintf('triple neg: %.2f, her2: %.2f, luminal: %.2f\n',100*mean(BM_at_roi(Dx==2)) ,100*mean(BM_at_roi(Dx==3)), 100*mean(BM_at_roi(Dx==1))) ;
    end
end
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(Pval);
[find(Pval<0.05) Pval(Pval<0.05) adj_p(Pval<0.05)]



vo_out = vo_atlas;
vo_out.fname='overlap_aal.nii';
vo_out.dt=[16,0];
IMG = zeros(vo_atlas.dim);
pctOverlap = 100*mean(CNT>unitVol);
for i=1:nroi,
   idx = find(atlas==i);
   IMG(idx)=pctOverlap(i);
end
spm_write_vol(vo_out,IMG);


return


% Crosstab analysis
%--------------------------------------------------------------------------
luminal = ordinal(T.luminal,{'none','luminal'},[],[0,0.5,1.5]);
HER2 = ordinal(T.HER2,{'none','HER2'},[],[0,0.5,1.5]);
TripleNega = ordinal(T.Basal,{'none','TripleNega'},[],[0,0.5,1.5]);



ODDs = zeros(nroi,1);
PVAL = zeros(nroi,1);
for i=1:nroi,
    X = [CNT(:,i)>unitVol];
    % [B,dev,stats] = mnrfit(X,luminal,'model','ordinal','interactions','off');
    % [B,dev,stats] = mnrfit(X,HER2,'model','ordinal','interactions','off');
    [B,dev,stats] = mnrfit(X,TripleNega,'model','ordinal','interactions','off');
    ODDs(i) = exp(B(2));
    PVAL(i) = stats.p(2);
end
[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(PVAL);
[find(PVAL<0.05) PVAL(PVAL<0.05) adj_p(PVAL<0.05)]


