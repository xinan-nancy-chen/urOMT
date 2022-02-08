function im_out = getData(data_tag,tind)
%this function loads the data images with no processing or interpolation
%last modified: 01/17/2017
%               01/23/2017
%               01/26/2017

addpath('../../Visualization/MATLAB_MRI_display/NIfTI_analyze/');

if nargin < 1
    data_tag = 'gauss';
    tind = 1;
elseif nargin < 2
    tind = 1;
end

cfg = getParams(data_tag);
if isfield(cfg,'dataset_name')
    switch cfg.dataset_name
        case {'MCA'} % 2D image
            im = double(imread(sprintf('%s%s%d%s',cfg.data_dir,cfg.data_name,tind,cfg.extension)));
            
        case {'LN','SD_LNC','ETOH','CohenAD','ATPKO','HN','test_dec21','microdia'} % 4D data
            D=load_untouch_nii(sprintf('%s%s%s',cfg.data_dir,cfg.data_name,cfg.extension));
            im = squeeze(D.img(:,:,:,tind));
        case {'people'}  % load mat
            im = load(sprintf('%s%s%d%s',cfg.data_dir,cfg.data_name,tind,cfg.extension));
            im = im.img;
            
        case {'Gauss'}
            a = 6; mv = 0.8;
            Nx = 20; Ny = Nx; Nz = Nx;
            x = linspace(-a,a,Ny) - 2*mv;
            y = linspace(-a,a,Nx) - 2*mv;
            z = linspace(-a,a,Nz) - 2*mv;
            [X,Y,Z] = meshgrid(x,y,z);
            rho0 = (100/sqrt(2*pi).*exp(-(X.^2/2)-(Y.^2/2)-(Z.^2/2)));
            x = x + mv; y = y + mv; z = z + mv;
            [X,Y,Z] = meshgrid(x,y,z);
            rho1 = (100/sqrt(2*pi).*exp(-(X.^2/2)-(Y.^2/2)-(Z.^2/2)));
            rho1 = imgaussfilt3(rho1,sqrt(0.2)*5) + rho1*0.1;
            rho.vol(1).data = rho0;
            rho.vol(2).data = rho1;
            im = rho.vol(tind).data;
            
        otherwise
            D=load_untouch_nii(sprintf('%s%s%02d%s%s',cfg.data_dir,cfg.data_name,tind,cfg.data_name_post,cfg.extension));
            im=D.img;
    end
            
else
    switch cfg.data_tag
        case {'tumor'}
            D=load_untouch_nii(sprintf('%s%s',cfg.data_dir,cfg.extension));
            im = squeeze(D.img(:,:,:,tind));
        otherwise
            D=load_untouch_nii(sprintf('%s%s%02d%s%s',cfg.data_dir,cfg.data_name,tind,cfg.data_name_post,cfg.extension));
            im=D.img;
    end
end

if length(size(im))==2
    im=double(im(cfg.x_range,cfg.y_range));
else
    im=double(im(cfg.x_range,cfg.y_range,cfg.z_range));
end

im(im<0) = 0;
im(im>10000) = 10000;

if isfield(cfg,'scalerho') && cfg.scalerho~=0
    im = im/cfg.scalerho;
end

if cfg.do_resize
    S_in = size(im);
    S_out = round(S_in*cfg.size_factor);
    im = resizeMatrix(im,S_out,'linear');%resizeMatrix(im,S_out,'cubic');
end


im_out = im;

end
