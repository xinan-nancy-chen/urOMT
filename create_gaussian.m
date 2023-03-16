% This script gives the creation of the gaussian spheres used in the
% example, as well as its chi functions.
addpath('utilities', 'Sensitivities','Inverse');

do_save = 1; % 1 if save the data to .mat

%% ====== 'gauss1' ======
tag = 'gauss1';
cfg = getParams(tag);
if cfg.mask_number == 0 % no mask is used in this case
    cfg.msk = ones(cfg.true_size);
end

% create gaussian sphere
a = 6; mv = 0.8;
Nx = 50; Ny = Nx; Nz = Nx;

x = linspace(-a,a,Ny) - 2*mv;
y = linspace(-a,a,Nx) - 2*mv;
z = linspace(-a,a,Nz) - 2*mv;

r_scale = [1,1.1,1.2,1.1,1];
chi_radius = 1.5;

for i = 1:(cfg.last_time-cfg.first_time)/cfg.time_jump+2
    [X,Y,Z] = meshgrid(x,y,z);
    tmp = (100/sqrt(2*pi).*exp(-(X.^2/2)-(Y.^2/2)-(Z.^2/2)));
    
    tmpchi = (X.^2/2)+(Y.^2/2)+(Z.^2/2)<chi_radius; 
    tmp(tmpchi>0) = r_scale(i)*tmp(tmpchi>0);
    cfg.tmpchi(i).data = tmpchi;
    figure, montageArray(tmpchi), axis image, clim([0,1]), title(sprintf('tmpchi%d',i))

    if cfg.do_resize
       tmp = resizeMatrix(tmp,round(cfg.size_factor.*size(tmp)),'linear');
    end
    if cfg.smooth>0
    tmp = affine_diffusion_3d(tmp,cfg.smooth,0.1,1,1);
    end
    
    % add diffusion
    if i~=1
    tmp = imgaussfilt3(tmp,sqrt(0.2)*i);
    end
    
    figure, montageArray(tmp), axis image, clim([0,50]), title(sprintf('rho%d',i))
    tmp(~cfg.msk) = 0;
    cfg.vol(i).data = tmp;
    x = x + mv; y = y + mv; z = z + mv;
    fprintf('i = %d, sum = %.5f\n', i, sum(tmp(:)))
    
end
%
% create chi function
x = linspace(-a,a,Ny) - 2*mv;
y = linspace(-a,a,Nx) - 2*mv;
z = linspace(-a,a,Nz) - 2*mv;

for i = 1:(cfg.last_time-cfg.first_time)/cfg.time_jump+1
    chi = zeros(prod(cfg.true_size),cfg.nt);
    for j = 1:cfg.nt
        [X,Y,Z] = meshgrid(x,y,z);
        tmpchi = (X.^2/2)+(Y.^2/2)+(Z.^2/2)<chi_radius; 
        chi(:,j) = tmpchi(:);
        x = x + mv/cfg.nt; y = y + mv/cfg.nt; z = z + mv/cfg.nt;
    end
    cfg.chi(i).data = chi;
    %figure, montageArray(chi), axis image, clim([0,1]), title(sprintf('chi%d',i))
end
%
if do_save
    save_dir_pre = sprintf('./data/%s/%s',cfg.dataset_name,tag);

    if ~exist(save_dir_pre,'dir')
        mkdir(save_dir_pre)
    end
    for i = 1:(cfg.last_time-cfg.first_time)/cfg.time_jump+2
        rho_n = cfg.vol(i).data;
        save(sprintf('%s/rho%d.mat',save_dir_pre,cfg.first_time+(i-1)*cfg.time_jump),'rho_n');
    end
    for i = 1:(cfg.last_time-cfg.first_time)/cfg.time_jump+1
        rho_n = cfg.chi(i).data;
        save(sprintf('%s/chi%d.mat',save_dir_pre,cfg.first_time+(i-1)*cfg.time_jump),'rho_n');
    end

end

%%







