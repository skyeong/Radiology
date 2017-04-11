function idmask = generate_rois(mnipos, RADIUS, fn_ref, fn_out)

% Load volume information
vref = spm_vol(fn_ref);

% Generate ROIs
IMG  = zeros(vref.dim);

% Generate Cubic Box
xx = mnipos(1)-RADIUS:mnipos(1)+RADIUS;
yy = mnipos(2)-RADIUS:mnipos(2)+RADIUS;
zz = mnipos(3)-RADIUS:mnipos(3)+RADIUS;
[X,Y,Z] = meshgrid(xx, yy, zz);
XYZ = [X(:), Y(:), Z(:), ones(size(Z(:)))];

% Converting MNI to Voxel space
Vxyz  = round(pinv(vref.mat)*XYZ');
try
    idroi = sub2ind(vref.dim, Vxyz(1,:), Vxyz(2,:), Vxyz(3,:));
catch
    idx = find(Vxyz(1,:)<1 | Vxyz(1,:)>vref.dim(1));
    idy = find(Vxyz(2,:)<1 | Vxyz(2,:)>vref.dim(2));
    idz = find(Vxyz(3,:)<1 | Vxyz(3,:)>vref.dim(3));
    out_of_FoV = union(idx,idy);
    out_of_FoV = union(out_of_FoV,idz);
    Vxyz(:,out_of_FoV) = [];
    idroi = sub2ind(vref.dim, Vxyz(1,:), Vxyz(2,:), Vxyz(3,:));
end

% Rounding Cubic-shape ROI to Sphere
for j=1:length(idroi),
    d = mnipos(1:3) - XYZ(j,1:3);
    radius = sqrt( sum(d.*d) );
    if radius <= RADIUS,
        IMG(idroi(j)) = 1;
    end
end
idmask = find(IMG>0);
% fprintf('center mni=[%d, %d, %d], nvox=%d)\n',mnipos,length(idmask));



% Save ROI-MASK image
if nargin>=4,
    vol = vref;
    vol.fname = fn_out;
    vol.dt = [4 0];
    vol.descrip = 'ROI';
    spm_write_vol(vol, IMG);
end