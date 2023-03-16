function im_out = getData(data_tag,data_tind,chi_tind)
%this function loads the data images and chi images with no processing or interpolation
%last modified: 03/04/2023 by Xinan Chen

%addpath('../../Visualization/MATLAB_MRI_display/NIfTI_analyze/');

if nargin < 1
    data_tag = 'gauss';
    data_tind = 1;
elseif nargin < 2
    data_tind = 1;
end

cfg = getParams(data_tag);
if isfield(cfg,'dataset_name')
    switch cfg.dataset_name
        case {'Gauss','RatBrainsCAA3M'}
            % load in data
            im_data = zeros([cfg.domain_size,length(data_tind)]);
            for j = 1:length(data_tind)
                im_data(:,:,:,j) = load_rhon(sprintf('%s%s%d%s',cfg.data_dir,cfg.data_name,data_tind(j),cfg.data_extension));
            end
            % load in chi if required
            if nargin > 2
                im_chi = zeros([cfg.domain_size,cfg.nt,length(chi_tind)]);
                for j = 1:length(chi_tind)
                    im_chi(:,:,:,:,j) = reshape(load_rhon(sprintf('%s%s%d%s',cfg.chi_dir,cfg.chi_name,chi_tind(j),cfg.chi_extension)),[cfg.domain_size,cfg.nt]);
                end
            end



        otherwise
            error('Unrecognized dataset name in getData.m!')
    end
            
else
    error('No dataset name is specified in getParams.m!')
end

if length(size(im_data))==2
    im_data=double(im_data(cfg.x_range,cfg.y_range));
elseif length(size(im_data))==3
    im_data=double(im_data(cfg.x_range,cfg.y_range,cfg.z_range));
else
    im_data=double(im_data(cfg.x_range,cfg.y_range,cfg.z_range,:));
end

im_data(im_data<0) = 0;
im_data(im_data>10000) = 10000;

if nargin > 2
if length(size(im_chi))==3
    im_chi=double(im_chi(cfg.x_range,cfg.y_range,:));
elseif length(size(im_chi))==4
    im_chi=double(im_chi(cfg.x_range,cfg.y_range,cfg.z_range,:));
else
    im_chi=double(im_chi(cfg.x_range,cfg.y_range,cfg.z_range,:,:));
end
end

im_out.im_data = im_data;
if nargin > 2
im_out.im_chi = im_chi;
end

end
