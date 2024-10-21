function cfg = getParams(data_name)
% This function returns parameters corresponding to the input 'data_name'

if nargin<1
    data_name='gaussian';
end
cfg.out_dir_pre = '../urOMT';

switch data_name
    case 'BC46_01' % for test
        cfg.data_tag = data_name;
        % input images
        cfg.data_dir = './data/BreastTumor/BC46_01/';
        cfg.data_name = 'BC46_01_frame';
        cfg.data_extension = '.mat';
        % input chi, comment if do not want to constrain source
        %cfg.chi_dir = './data/Gauss/gauss3/';
        %cfg.chi_name = 'chi';
        %cfg.chi_extension = '.mat';

        %cfg.anato = './data/BreastTumor/BC46_01/BC46_01_anatomy.mat';

        cfg.dataset_name = 'BreastTumor';

        cfg.mask_number = 1;
        switch cfg.mask_number
            case 1 % tumor ROI mask (dilated by 2)
                cfg.ROI_msk_path = './data/BreastTumor/BC46_01/mask4.mat';
                cfg.x_range = 80:112;
                cfg.y_range = 223:253;
                cfg.z_range = 49:77;
        end
        
        cfg.domain_size = [320,320,128];
        cfg.spacing.x = 1.0625; 
        cfg.spacing.y = 1.0625; 
        cfg.spacing.z = 1.4;

        cfg.n1=length(cfg.x_range);
        cfg.n2=length(cfg.y_range);
        cfg.n3=length(cfg.z_range);
        
        cfg.do_resize = 0;
        cfg.size_factor = 1;
        
        cfg.true_size=round(cfg.size_factor*[length(cfg.x_range),length(cfg.y_range),length(cfg.z_range)]);
        
        cfg.reinitR = 0;%0 if do consecutively and 1 if reinitialize rho
        cfg.smooth = 1;
        
        cfg.first_time = 5;
        cfg.time_jump = 2;
        cfg.last_time = 25;
        
        cfg.sigma = 2e-3;
        cfg.sig_str = '2e3';
        cfg.dt = 0.4;
        cfg.nt = 10;
        cfg.alpha = 30000; % weight for source
        cfg.beta  = 1000; % weight for matching final image
        cfg.niter_pcg = 60;
        
        cfg.dTri = 1; %1 := 'closed', 3:= 'open'

        cfg.version = sprintf('diff_%s_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_beta_%d_alpha_%d_smooth%d_dtri%d_rreinit%d_pcg%d',...
            cfg.sig_str,cfg.time_jump,cfg.dt,cfg.nt,cfg.first_time,cfg.last_time,cfg.beta,cfg.alpha,cfg.smooth,cfg.dTri,cfg.reinitR,cfg.niter_pcg);
        cfg.out_dir=sprintf('%s/%s/%s/%s',cfg.out_dir_pre,cfg.dataset_name,cfg.data_tag,cfg.version);
        %%

    case 'C1217' % for test
        cfg.data_tag = data_name;
        % input images
        cfg.data_dir = './data/RatBrainsCAA3M/C1217/';
        cfg.data_name = 'C1217_psnrv';
        cfg.data_extension = '.mat';
        % input chi, comment if do not want to constrain source
        %cfg.chi_dir = './data/Gauss/gauss3/';
        %cfg.chi_name = 'chi';
        %cfg.chi_extension = '.mat';

        cfg.anato = './data/RatBrainsCAA3M/C1217/C1217_anatomy.mat';

        cfg.dataset_name = 'RatBrainsCAA3M';

        cfg.mask_number = 1;
        switch cfg.mask_number
            case 1 %brain-only mask
                cfg.ROI_msk_path = './data/RatBrainsCAA3M/C1217/mask1.mat';
                cfg.x_range = 24:79;
                cfg.y_range = 1:106;
                cfg.z_range = 31:81;
        end
        
        cfg.domain_size = [100,106,100];
        cfg.spacing.x = 1; 
        cfg.spacing.y = 1; 
        cfg.spacing.z = 1;

        cfg.n1=length(cfg.x_range);
        cfg.n2=length(cfg.y_range);
        cfg.n3=length(cfg.z_range);
        
        cfg.do_resize = 0;
        cfg.size_factor = 1;
        
        cfg.true_size=round(cfg.size_factor*[length(cfg.x_range),length(cfg.y_range),length(cfg.z_range)]);
        
        cfg.reinitR = 0;%0 if do consecutively and 1 if reinitialize rho
        cfg.smooth = 1;
        
        cfg.first_time = 7;
        cfg.time_jump = 2;
        cfg.last_time = 33;
        
        cfg.sigma = 2e-3;
        cfg.sig_str = '2e3';
        cfg.dt = 0.4;
        cfg.nt = 10;
        cfg.alpha = 10000; % weight for source
        cfg.beta  = 50; % weight for matching final image
        cfg.niter_pcg = 60;
        
        cfg.dTri = 1; %1 := 'closed', 3:= 'open'

        cfg.version = sprintf('diff_%s_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_beta_%d_alpha_%d_smooth%d_dtri%d_rreinit%d_pcg%d',...
            cfg.sig_str,cfg.time_jump,cfg.dt,cfg.nt,cfg.first_time,cfg.last_time,cfg.beta,cfg.alpha,cfg.smooth,cfg.dTri,cfg.reinitR,cfg.niter_pcg);
        cfg.out_dir=sprintf('%s/%s/%s/%s',cfg.out_dir_pre,cfg.dataset_name,cfg.data_tag,cfg.version);
        %%
    case 'gauss1' % for test
        cfg.data_tag = data_name;
        % input images
        cfg.data_dir = './data/Gauss/gauss1/';
        cfg.data_name = 'rho';
        cfg.data_extension = '.mat';
        % input chi, comment if do not want to constrain source
        cfg.chi_dir = './data/Gauss/gauss1/';
        cfg.chi_name = 'chi';
        cfg.chi_extension = '.mat';

        cfg.dataset_name = 'Gauss';

        cfg.mask_number = 0;
        switch cfg.mask_number
            case 0 % use no mask
                cfg.x_range = 1:50;
                cfg.y_range = 1:50;
                cfg.z_range = 1:50;
            case 1 % if use any mask of region of interest
                cfg.ROI_msk_path = '';
                cfg.x_range = 1:50;
                cfg.y_range = 1:50;
                cfg.z_range = 1:50;
        end
        
        cfg.domain_size = [50,50,50];
        cfg.spacing.x = 1; 
        cfg.spacing.y = 1; 
        cfg.spacing.z = 1;

        cfg.n1=length(cfg.x_range);
        cfg.n2=length(cfg.y_range);
        cfg.n3=length(cfg.z_range);
        
        cfg.do_resize = 0;
        cfg.size_factor = 1;
        
        cfg.true_size=round(cfg.size_factor*[length(cfg.x_range),length(cfg.y_range),length(cfg.z_range)]);
        
        cfg.reinitR = 1;%0 if do consecutively and 1 if reinitialize rho
        cfg.smooth = 0;
        
        cfg.first_time = 1;
        cfg.time_jump = 1;
        cfg.last_time = 4;
        
        cfg.sigma = 2e-3;
        cfg.sig_str = '2e3';
        cfg.dt = 0.4;
        cfg.nt = 10;
        cfg.alpha = 9000; % weight for source
        cfg.beta  = 5000; % weight for matching final image
        cfg.niter_pcg = 60;
        
        cfg.dTri = 1; %1 := 'closed', 3:= 'open'

        cfg.version = sprintf('diff_%s_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_beta_%d_alpha_%d_smooth%d_dtri%d_rreinit%d_pcg%d',...
            cfg.sig_str,cfg.time_jump,cfg.dt,cfg.nt,cfg.first_time,cfg.last_time,cfg.beta,cfg.alpha,cfg.smooth,cfg.dTri,cfg.reinitR,cfg.niter_pcg);
        cfg.out_dir=sprintf('%s/%s/%s/%s',cfg.out_dir_pre,cfg.dataset_name,cfg.data_tag,cfg.version);
        %%
    otherwise
        disp('data-name unrecognized')
end
