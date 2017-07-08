vref = spm_vol('/Users/skyeong/connectome/atlas/faal_s.nii');
ROI = spm_read_vols(vref);

fn_out = 'AAL';
pctOverlap = round(pctOverlap);
mylist = unique(pctOverlap(:));
nCi = length(mylist);

IMG = zeros(vref.dim);
IMG_L = zeros(vref.dim);
IMG_R = zeros(vref.dim);
for i=1:nCi,
    idroi = find(abs(pctOverlap-mylist(i))<0.1);
    idx = [];
    idx_L = [];
    idx_R = [];
    for j=1:length(idroi),
        idx = [idx; find(ROI==idroi(j))];
        if idroi(j)>=109,
            idx_L = [idx_L; find(ROI==idroi(j))];
        elseif rem(idroi(j),2)==1,  % for odd (left)
            idx_L = [idx_L; find(ROI==idroi(j))];
        else
            idx_R = [idx_R; find(ROI==idroi(j))];
        end
    end
    IMG(idx)=mylist(i);
    IMG_L(idx_L)=mylist(i);
    IMG_R(idx_R)=mylist(i);
end


vout = vref;
vout.dt = [16 0];
vout.fname = [fn_out '.nii'];
spm_write_vol(vout,IMG);

vout = vref;
vout.dt = [16 0];
vout.fname = [fn_out '_L.nii'];
spm_write_vol(vout,IMG_L);

vout = vref;
vout.dt = [16 0];
vout.fname = [fn_out '_R.nii'];
spm_write_vol(vout,IMG_R);