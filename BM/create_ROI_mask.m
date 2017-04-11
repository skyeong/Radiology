function create_ROI_mask(MNIcenter, radius, fn_out, fn_mask,ROIshape)

if nargin<5,
    ROIshape='sphere';
end


%  GENERATE SEED ROIS DEFINED ABOVE
%__________________________________________________________________________

vref = spm_vol(fn_mask);

Vxyz=[];
for i=-radius:radius,
    for j=-radius:radius,
        for k=-radius:radius,
            RADIUS = sqrt(i*i + j*j + k*k);
            if strcmpi(ROIshape,'sphere'),
                if (RADIUS <= radius)
                    sROI1 = [MNIcenter(1)+i,MNIcenter(2)+j,MNIcenter(3)+k];
                    Vxyz = [Vxyz;sROI1];
                end
            else
                sROI1 = [MNIcenter(1)+i,MNIcenter(2)+j,MNIcenter(3)+k];
                Vxyz = [Vxyz;sROI1];
            end
        end
    end
end
Vxyz = [Vxyz, ones(size(Vxyz,1),1)];
Rxyz = round(inv(vref.mat)*Vxyz')';
idroi = sub2ind(vref.dim(1:3),Rxyz(:,1),Rxyz(:,2),Rxyz(:,3));
idroi = unique(idroi);


%  SAVE ROI AS NIFTI FILE FORMAT
%__________________________________________________________________________

vo = vref;
vo.fname = fn_out;
IMG = zeros(vref.dim);
IMG(idroi) = radius;
spm_write_vol(vo,IMG);

fprintf('ROI is saved at\n');
fprintf('%s\n',fn_out);


