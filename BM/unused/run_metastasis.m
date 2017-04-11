warning('off','all');
addpath('/Users/skyeong/matlabwork/spm12');


Ithr = 650;
DATApath = '/Volumes/JetDrive/data/SJ_Ahn';
subjlist = {'0707021232','0707021510','133381','6887227','6898853','806160','6896997','8273822'};
nsubj = length(subjlist);

% Load CSF mask image;
fn_mask = fullfile(DATApath,'brainmask','rCSF.nii');
vo_mask = spm_vol(fn_mask);
MASK = spm_read_vols(vo_mask);
idmask = find(MASK>0);

% Extract abnormal tissues
for c=1:nsubj
% for c=2
    subjname = subjlist{c};
    
    % Read metastasis image
    fn = fullfile(DATApath,'metastasis',subjname,'wmetastasis_bet.nii');
    vo = spm_vol(fn);
    M = spm_read_vols(vo);
    DIM = vo.dim;
    
    % Create Empty Images
    SIG = zeros(DIM);
    BKG = zeros(DIM);
    
    % Find indices with hyperintensities (>threshold)
    idx = find(M>Ithr);
    

    % Spliting into clusters
    [vx,vy,vz] = ind2sub(DIM, idx);
    Vxyz = [vx vy vz ones(size(vx))];
    ilab = spm_clusters(Vxyz(:,1:3)');
    nCl = max(ilab);
    
    % Identifying cancers and noises
    S = struct();
    B = struct();
    cnt_s=1;
    cnt_b=1;
    for i=1:nCl,
        idx_cl = find(ilab==i);
        if length(idx_cl)<3, continue; end
        Rxyz = vo.mat*Vxyz(idx_cl,:)';
        muRxyz = mean(Rxyz(1:3,:),2);
        
        
        % Calc properties
        [fa, fa2] = CalcProperties(Rxyz(1:3,:));
        
        
%         aa = muRxyz-[-21.2 -76.2 -14.0]';
%         aa = sqrt(sum(aa.*aa));
%         if aa<10,
%             fa
%             fa2
%         end
        
        % Compute overlap with CSF mask
        idx_ovp = intersect(idx(idx_cl),idmask);
        ovlpct = length(idx_ovp)./length(idx_cl);
        if ovlpct>=0.20,
            BKG(idx(idx_cl)) = 1;
            B(cnt_b).fa = fa;
            B(cnt_b).fa2 = fa2;
            B(cnt_b).xyz = muRxyz;
            cnt_b = cnt_b + 1;
            continue;
        end
        
    
        % Predefined noise components (around temporal poles)
        if isFalsePositive(muRxyz),
            BKG(idx(idx_cl)) = 1;
            B(cnt_b).fa = fa;
            B(cnt_b).fa2 = fa2;
            B(cnt_b).xyz = muRxyz;
            cnt_b = cnt_b + 1;
            continue;
        end
        
        if muRxyz(2)>-30 && length(idx_cl)>20,
            SIG(idx(idx_cl)) = 1;
            S(cnt_s).fa = fa;
            S(cnt_s).fa2 = fa2;
            S(cnt_s).xyz = muRxyz;
            cnt_s = cnt_s + 1;
        elseif fa<eps || fa2<eps,
            BKG(idx(idx_cl)) = 1;
            B(cnt_b).fa = fa;
            B(cnt_b).fa2 = fa2;
            B(cnt_b).xyz = muRxyz;
            cnt_b = cnt_b + 1;
        elseif abs(fa-fa2)<0.2,
            BKG(idx(idx_cl)) = 1;
            B(cnt_b).fa = fa;
            B(cnt_b).fa2 = fa2;
            B(cnt_b).xyz = muRxyz;
            cnt_b = cnt_b + 1;
        elseif fa<0.8 && fa2<0.7,
            SIG(idx(idx_cl)) = 1;
            S(cnt_s).fa = fa;
            S(cnt_s).fa2 = fa2;
            S(cnt_s).xyz = muRxyz;
            cnt_s = cnt_s + 1;
        else
            BKG(idx(idx_cl)) = 1;
            B(cnt_b).fa = fa;
            B(cnt_b).fa2 = fa2;
            B(cnt_b).xyz = muRxyz;
            cnt_b = cnt_b + 1;
        end
        
 
    end
    
    % Save resulting output
    vout = vo;
    vout.fname=fullfile(DATApath,'metastasis',subjname,'SIG.nii');
    spm_write_vol(vout,SIG);
    
    vout.fname=fullfile(DATApath,'metastasis',subjname,'BKG.nii');
    spm_write_vol(vout,BKG);
    
    % figure; hist(FA2(FA2>0),20);
    % figure; hist(FA(FA>0),20);
    % pause
end