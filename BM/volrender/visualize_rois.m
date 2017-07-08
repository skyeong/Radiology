function h0=visualize_rois(volname,idx,col,overlay,outname)
% visualize_rois(volname,idx,col,overlay,outname)
%
% Renders file specified by volname in 3D
% volname in string, idx specifies region index to plot as separate
% col is colormap, overlay for brain overlay (0 to turn off, 1 for
% left hemisphere, 2 for right hemisphere, otherwise both hemisphere),
% and outname for output figure file
%
% First written by Bumhee Park, 2011
% Functions compiled by Jeong Hoon Ko, 6/21/2012

if nargin<1; volname=spm_select(1,'any'); end
vol=spm_vol(volname);
origvol=spm_read_vols(vol);
if nargin<5;
    outname=[]; printflag=0;
else
    printflag=1;
    if ~strcmp(outname(end-3:end),'.png')
        outname=[outname '.png'];
    end
end
if nargin<4; overlay=0; end
if nargin<2; idx=1:max(origvol(:)); end
if nargin<3; col=jet(max(idx)); end
if length(idx)>size(col,1); error('Provided colormap is too small'); end

rdidx=idx;%randperm(length(idx));
%h0=figure;
%set(h0,'Position',[500 300 560 420]);
%h=waitbar(0,'Please wait...');
%set(h,'Position',[100 300 270 56.25]);
hold on
for i=1:length(idx)
    % waitbar(i/length(idx))
    %disp(num2str(i));
    now=zeros(size(origvol));
    now(origvol==idx(i))=1;
    if(length(find(origvol==idx(i)))>10)
        a=spm_surf_mod(vol.fname,2,0.2,now);
        %a=spm_surf_mod(vol.fname,2,1,now);
        s=gifti(a.surffile);
        s0=s;
        s1=spm_mesh_inflate(s0,2);
        fv.faces=s1.faces;
        fv.vertices=s1.vertices;
        opt.light=0;
        hsk=vol_rendobj_mod(fv,col(i,:),opt);
        %         set(hsk,'FaceAlpha',0.3);
    end
    %hsk=vol_rendobj_mod(fv,col(rdidx(i),:),opt);
    
end
%close(h)




if overlay
    %brain overlay
    s=gifti([spm('dir') '/canonical/cortex_5124.surf.gii']);
    %inflate brain overlay to match AAL; perform separate for Lt and Rt
    lv.vertices=s.vertices(1:round(size(s.vertices,1)/2),:);
    rv.vertices=s.vertices(round(size(s.vertices,1)/2)+1:end,:);
    lv=spheroid_inflate(lv,5);
    rv=spheroid_inflate(rv,5);
    fv.vertices=[lv.vertices; rv.vertices];
    fv.faces=s.faces;
    fv=spm_mesh_inflate(fv,3);
    switch overlay
        case 1
            fv.faces=fv.faces(1:round(size(fv.faces,1)/2),:); %only left hemisphere
        case 2
            fv.faces=fv.faces(round(size(fv.faces,1)/2)+1:end,:); %only right hemisphere
        otherwise
            fv.faces=fv.faces; % %whole brain
    end
    hsk=vol_rendobj_mod(fv,[1 1 1].*0.75); set(hsk,'FaceAlpha',0.3);
end

%post processing of rendering
material dull
lighting phong;
% camlight(90,-10);
% camlight(-90,10);

view(3); axis image;
view(-90,10);
axis off;
hold off

set(gcf,'InvertHardCopy','off','Color',[1 1 1]);
if printflag,
    print('-dpng','-r600',outname)
end

end

function s0=spheroid_inflate(s0,m)
%inflates s0.vertices in outward direction from the centroid
%by magnitude m
cent=mean(s0.vertices,1);
n=s0.vertices-repmat(cent,size(s0.vertices,1),1);
lenn=sum(n.^2,2).^0.5;
n=n./repmat(lenn,1,3);
s0.vertices=s0.vertices+m*n;
end

function varargout = spm_surf_mod(P,mode,thresh,ovol)
% Surface extraction
% FORMAT spm_surf(P,mode,thresh)
%
% P      - char array of filenames
%          Usually, this will be c1xxx.img & c2xxx.img - grey and white
%          matter segments created using the segmentation routine.
% mode   - operation mode [1: rendering, 2: surface, 3: both]
% thresh - vector or threshold values for extraction [default: 0.5]
%          This is only relevant for extracting surfaces, not rendering.
%
% Generated files (depending on 'mode'):
% A "render_xxx.mat" file can be produced that can be used for
% rendering activations on to, see spm_render.
%
% A "xxx.surf.gii" file can also be written, which is created using
% Matlab's isosurface function.
% This extracted brain surface can be viewed using code something like:
%    FV = gifti(spm_select(1,'mesh','Select surface data'));
%    FV = export(FV,'patch');
%    fg = spm_figure('GetWin','Graphics');
%    ax = axes('Parent',fg);
%    p  = patch(FV, 'Parent',ax,...
%           'FaceColor', [0.8 0.7 0.7], 'FaceVertexCData', [],...
%           'EdgeColor', 'none',...
%           'FaceLighting', 'phong',...
%           'SpecularStrength' ,0.7, 'AmbientStrength', 0.1,...
%           'DiffuseStrength', 0.7, 'SpecularExponent', 10);
%    set(0,'CurrentFigure',fg);
%    set(fg,'CurrentAxes',ax);
%    l  = camlight(-40, 20);
%    axis image;
%    rotate3d on;
%
% FORMAT out = spm_surf(job)
%
% Input
% A job structure with fields
% .data   - cell array of filenames
% .mode   - operation mode
% .thresh - thresholds for extraction
% Output
% A struct with fields (depending on operation mode)
% .rendfile - cellstring containing render filename
% .surffile - cellstring containing surface filename(s)
%__________________________________________________________________________
%
% This surface extraction is not particularly sophisticated.  It simply
% smooths the data slightly and extracts the surface at a threshold of
% 0.5. The input segmentation images can be manually cleaned up first using
% e.g., MRIcron.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_surf.m 4185 2011-02-01 18:46:18Z guillaume $

%

%SVNrev = '$Rev: 4185 $';

%spm('FnBanner',mfilename,SVNrev);
%spm('FigName','Surface');

%-Get input: filenames 'P'
%--------------------------------------------------------------------------
try
    if isstruct(P)
        job    = P;
        P      = strvcat(job.data);
        mode   = job.mode;
        thresh = job.thresh;
    end
catch
    [P, sts] = spm_select([1 Inf],'image','Select images');
    if ~sts, varargout = {}; return; end
end

%-Get input: operation mode 'mode'
%--------------------------------------------------------------------------
try
    mode;
catch
    mode = spm_input('Save','+1','m',...
        ['Save Rendering|'...
        'Save Extracted Surface|'...
        'Save Rendering and Surface'],[1 2 3],3);
end

%-Get input: threshold for extraction 'thresh'
%--------------------------------------------------------------------------
try
    thresh;
catch
    thresh = 0.5;
end

%-Surface extraction
%--------------------------------------------------------------------------
spm('FigName','Surface: working');
spm('Pointer','Watch');
try
    out = do_it(P,mode,thresh,ovol);
catch
    out = do_it(P,mode,thresh);
end
spm('Pointer','Arrow');
spm('FigName','Surface: done');

if nargout > 0
    varargout{1} = out;
end

return;
end
%==========================================================================
function out = do_it(P,mode,thresh,ovol)

V  = spm_vol(P);
br = zeros(V(1).dim(1:3));

if exist('ovol')
    br=ovol;
else
    for i=1:V(1).dim(3),
        B         = spm_matrix([0 0 i]);
        tmp       = spm_slice_vol(V(1),B,V(1).dim(1:2),1);
        for j=2:length(V),
            M   = V(j).mat\V(1).mat*B;
            tmp = tmp + spm_slice_vol(V(j),M,V(1).dim(1:2),1);
        end
        br(:,:,i) = tmp;
    end
end

% Build a 3x3x3 seperable smoothing kernel and smooth
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;
spm_conv_vol(br,br,kx,ky,kz,-[1 1 1]);

[pth,nam,ext] = fileparts(V(1).fname);

if any(mode==[1 3])
    % Produce rendering
    %----------------------------------------------------------------------
    out.rendfile{1} = fullfile(pth,['render_' nam '.mat']);
    tmp = struct('dat',br,'dim',size(br),'mat',V(1).mat);
    renviews(tmp,out.rendfile{1});
end

if any(mode==[2 3])
    % Produce extracted surface
    %----------------------------------------------------------------------
    for k=1:numel(thresh)
        [faces,vertices] = isosurface(br,thresh(k));
        
        % Swap around x and y because isosurface does for some
        % wierd and wonderful reason.
        Mat      = V(1).mat(1:3,:)*[0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];
        vertices = (Mat*[vertices' ; ones(1,size(vertices,1))])';
        if numel(thresh)==1
            nam1 = nam;
        else
            nam1 = sprintf('%s-%d',nam,k);
        end
        out.surffile{k} = fullfile(pth,[nam1 '.surf.gii']);
        save(gifti(struct('faces',faces,'vertices',vertices)),out.surffile{k});
    end
end
return;
end

%==========================================================================
function renviews(V,oname)
% Produce images for rendering activations to
%
% FORMAT renviews(V,oname)
% V     - mapped image to render, or alternatively
%         a structure of:
%         V.dat - 3D array
%         V.dim - size of 3D array
%         V.mat - affine mapping from voxels to millimeters
% oname - the name of the render.mat file.
%__________________________________________________________________________
%
% Produces a matrix file "render_xxx.mat" which contains everything that
% "spm_render" is likely to need.
%
% Ideally, the input image should contain values in the range of zero
% and one, and be smoothed slightly.  A threshold of 0.5 is used to
% distinguish brain from non-brain.
%__________________________________________________________________________

linfun = inline('fprintf([''%-30s%s''],x,[repmat(sprintf(''\b''),1,30)])','x');
linfun('Rendering: ');

linfun('Rendering: Transverse 1..'); rend{1} = make_struct(V,[pi 0 pi/2]);
linfun('Rendering: Transverse 2..'); rend{2} = make_struct(V,[0 0 pi/2]);
linfun('Rendering: Saggital 1..');   rend{3} = make_struct(V,[0 pi/2 pi]);
linfun('Rendering: Saggital 2..');   rend{4} = make_struct(V,[0 pi/2 0]);
linfun('Rendering: Coronal 1..');    rend{5} = make_struct(V,[pi/2 pi/2 0]);
linfun('Rendering: Coronal 2..');    rend{6} = make_struct(V,[pi/2 pi/2 pi]);

linfun('Rendering: Save..');
if spm_check_version('matlab','7') >= 0
    save(oname,'-V6','rend');
else
    save(oname,'rend');
end
linfun('                 ');
if ~spm('CmdLine')
    disp_renderings(rend);
    spm_print;
end
return;
end

%==========================================================================
function str = make_struct(V,thetas)
[D,M]     = matdim(V.dim(1:3),V.mat,thetas);
[ren,dep] = make_pic(V,M*V.mat,D);
str       = struct('M',M,'ren',ren,'dep',dep);
return;
end

%==========================================================================
function [ren,zbuf] = make_pic(V,M,D)
% A bit of a hack to try and make spm_render_vol produce some slightly
% prettier output.  It kind of works...
if isfield(V,'dat'), vv = V.dat; else vv = V; end;
[REN, zbuf, X, Y, Z] = spm_render_vol(vv, M, D, [0.5 1]);
fw        = max(sqrt(sum(M(1:3,1:3).^2)));
msk       = find(zbuf==1024);
brn       = ones(size(X));
brn(msk)  = 0;
brn       = spm_conv(brn,fw);
X(msk)    = 0;
Y(msk)    = 0;
Z(msk)    = 0;
msk       = find(brn<0.5);
tmp       = brn;
tmp(msk)  = 100000;
sX        = spm_conv(X,fw)./tmp;
sY        = spm_conv(Y,fw)./tmp;
sZ        = spm_conv(Z,fw)./tmp;
zbuf      = spm_conv(zbuf,fw)./tmp;
zbuf(msk) = 1024;

vec       = [-1 1 3]; % The direction of the lighting.
vec       = vec/norm(vec);
[t,dx,dy,dz] = spm_sample_vol(vv,sX,sY,sZ,3);
IM        = inv(diag([0.5 0.5 1])*M(1:3,1:3))';
ren       = IM(1:3,1:3)*[dx(:)' ; dy(:)' ; dz(:)'];
len       = sqrt(sum(ren.^2,1))+eps;
ren       = [ren(1,:)./len ; ren(2,:)./len ; ren(3,:)./len];
ren       = reshape(vec*ren,[size(dx) 1]);
ren(ren<0) = 0;
ren(msk)  = ren(msk)-0.2;
ren       = ren*0.8+0.2;
mx        = max(ren(:));
ren       = ren/mx;
return;
end

%==========================================================================
function Fgraph = disp_renderings(rend)
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);
hght = 0.95;
nrow = ceil(length(rend)/2);
ax=axes('Parent',Fgraph,'units','normalized','Position',[0, 0, 1, hght],'Visible','off');
image(0,'Parent',ax);
set(ax,'YTick',[],'XTick',[]);

for i=1:length(rend),
    ren = rend{i}.ren;
    ax=axes('Parent',Fgraph,'units','normalized',...
        'Position',[rem(i-1,2)*0.5, floor((i-1)/2)*hght/nrow, 0.5, hght/nrow],...
        'Visible','off');
    image(ren*64,'Parent',ax);
    set(ax,'DataAspectRatio',[1 1 1], ...
        'PlotBoxAspectRatioMode','auto',...
        'YTick',[],'XTick',[],'XDir','normal','YDir','normal');
end
drawnow;
return;
end

%==========================================================================
function [d,M] = matdim(dim,mat,thetas)
R = spm_matrix([0 0 0 thetas]);
bb = [[1 1 1];dim(1:3)];
c  = [  bb(1,1) bb(1,2) bb(1,3) 1
    bb(1,1) bb(1,2) bb(2,3) 1
    bb(1,1) bb(2,2) bb(1,3) 1
    bb(1,1) bb(2,2) bb(2,3) 1
    bb(2,1) bb(1,2) bb(1,3) 1
    bb(2,1) bb(1,2) bb(2,3) 1
    bb(2,1) bb(2,2) bb(1,3) 1
    bb(2,1) bb(2,2) bb(2,3) 1]';
tc = diag([2 2 1 1])*R*mat*c;
tc = tc(1:3,:)';
mx = max(tc);
mn = min(tc);
M  = spm_matrix(-mn(1:2))*diag([2 2 1 1])*R;
d  = ceil(abs(mx(1:2)-mn(1:2)))+1;
return;
end

function [hiso,FV,N]=vol_rendobj_mod(vol,isoval,opt,fcol,mode)
% on Feb., 23, 2003, by Hae-Jeong Park
% modified by Jeong Hoon Ko
% 6/21/2012: added as function in visualize_rois, added opt.light
% to remove standard lighting option for customization in the main fcn

hiso=[]; FV=[]; N=[];

if nargin<1,
    vol=vol_vol;
end;
if nargin<2, isoval=[]; end;
if nargin<3, opt.light=1; end;
if nargin<4, fcol=[0 0.5 0.5]; end;
if nargin<5, mode='nobar'; end;


R=[];
if ischar(vol),
    v=spm_vol(vol);
    R=spm_read_vols(v);
    id=find(R>0);
    isoval=mean(R(id));
    id=find(R>=isoval);
    isoval=0.9*mean(R(id));
    FV=isosurface(R,isoval);
elseif isnumeric(vol),
    R=vol; id=find(R>0);
    if isempty(isoval),
        isoval=mean(R(id));
        id=find(R>=isoval);
        isoval=0.9*mean(R(id));
    end;
    FV=isosurface(R,isoval);
elseif isstruct(vol),
    if isfield(vol,'vertices'), %vol_rendobj(FV,fcol,mode)
        FV=vol;
        if nargin>=2, fcol=isoval; end;
    else
        R=spm_read_vols(vol);
        id=find(R>0);
        isoval=mean(R(id));
        id=find(R>=isoval);
        isoval=0.9*mean(R(id));
        FV=isosurface(R,isoval);
    end;
end;
%    FV.vertices=FV.vertices(:,[2 1 3]);
%    fprintf('vertex:%d\n',size(FV.vertices,1));

a=max(size(FV.vertices));
b=max(size(fcol));
N=[];
if isfield(FV,'normals'),
    N=FV.normals;
    FV=rmfield(FV,'normals');
end;
if a~=b,
    hiso = patch(FV,'FaceColor', fcol, 'EdgeColor', 'none');
else,
    if min(size(fcol))==3,
        if max(fcol(:))>1,fcol=fcol/255; end;
        hiso = patch(FV,  'FaceVertexCData',fcol,'EdgeColor', 'none','CDataMapping','direct','FaceColor', 'interp');
    else
        colormap(jet(128));
        fcol=fcol(:);
        hiso = patch(FV,  'FaceVertexCData',fcol, 'FaceColor', 'interp','EdgeColor', 'none');
    end;
end;

if isempty(R)==0,
    FV.normal=isonormals(R,hiso);
end;
if isempty(N),
else
    FV.normals=N;
end;

if isfield(FV,'normals'),   set(hiso, 'vertexnormals', FV.normals); set(hiso,'AmbientStrength',.8);end;
% hcap = patch(isocaps(vol,isoval), 'FaceColor',fcol, 'EdgeColor', 'none'); %'interp'
% set(hcap,'AmbientStrength',.8);

set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50);
if opt.light,
    view(3); axis image;
    camlight(90,-10); camlight(-90,10);
    lighting phong
    material dull;
    axis off;
    view(180,90);
end
if strcmp(mode,'bar'),
    h1=colorbar('vert');
    pos=get(h1,'Position');
    pos(2)=pos(2)*1.5;
    pos(4)=pos(4)*0.8;
    set(h1,'Position',pos);
end;
end