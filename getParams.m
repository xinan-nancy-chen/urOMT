function cfg = getParams(data_name)

% This function returns parameters corresponding to the input 'data_name'

if nargin<1
    data_name='gaussian';
end
cfg.out_dir_pre = '../urOMT';

switch data_name
    %% 3-month test
    case '022619A_psnrv'%%C1217 (WT)
        cfg.data_tag = '022619A_psnrv';
        cfg.data_dir = '../../data/3months_data/psnrv/WT/C1217/psnrv/';
        cfg.data_name = 'psnrv_C1217_022619A_DOTA37_30ul_E';
        cfg.data_name_post = '';
        cfg.extension = '.nii';
        
        cfg.max_dpsnrv = '../../data/3months_data/MAXpsnrv/C1217_022619A_psnrv_max.nii';
        cfg.anato = '../../data/3months_data/psnrv/WT/C1217/pbase_snrv_C1217_022619A_DOTA37_30ul_E48.nii';
        
        cfg.sp_mask_opts(1).name = 'brain'; % brain tissue = 2, CSF = 3, take mask as 2 and 3
        cfg.sp_mask_opts(1).path = '../../data/3months_data/3months_mask_brainCSF/C1217.nii';
        cfg.sp_mask_opts(1).type(1).label = 'CSF';
        cfg.sp_mask_opts(1).type(1).value = 3;
        cfg.sp_mask_opts(1).type(2).label = 'tissue';
        cfg.sp_mask_opts(1).type(2).value = 2;
        
        cfg.dataset_name = 'CAA';
        cfg.name = 'C1217';
        cfg.mask_number = 1;
        switch cfg.mask_number
            case 1%head mask
                cfg.ROI_msk_path = '../../data/3months_data/maskforOMT/C1217_mask.nii';
                cfg.x_range = 24:79;%26:77;
                cfg.y_range = 1:106;%1:106;
                cfg.z_range = 31:81;%33:79;
        end
        
        cfg.domain_size = [100,106,100];
        
        cfg.n1=length(cfg.x_range);
        cfg.n2=length(cfg.y_range);
        cfg.n3=length(cfg.z_range);
        
        cfg.do_resize = 0;
        cfg.size_factor = 1;
        
        cfg.scalerho = 0;%100;%0;
        if cfg.scalerho==0
           scale_str = '';
        else
           scale_str = sprintf('_scale%d',cfg.scalerho);
        end
        
        cfg.data_index_E = 14:48;

        cfg.true_size=round(cfg.size_factor*[length(cfg.x_range),length(cfg.y_range),length(cfg.z_range)]);
        
        cfg.smooth = 1;
        cfg.reinitR = 0;%0 if do consecutively and 1 if reinitialize rho
        
        cfg.first_time = cfg.data_index_E(7);%cfg.data_index_E(5);%cfg.data_index_E(13);
        cfg.time_jump = 2;
        cfg.last_time = cfg.data_index_E(33);%cfg.data_index_E(27);%cfg.data_index_E(25);%cfg.data_index_E(33);
        
        cfg.sigma = 2e-3;
        cfg.sig_str = '2e3';
        cfg.dt = 0.4;
        cfg.nt = 10;
        cfg.alpha = 3000;%10000;%50000;%8000;%1000;%100;%1000; % weight for source
        cfg.beta  = 50;%5000;%500;%5000; % weight for matching final image
        cfg.niter_pcg = 60;
        
        cfg.dTri = 3;%1 := 'closed', 3:= 'open'
        
        cfg.version = sprintf('diff_%s_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_beta_%d_alpha_%d_smooth%d_dtri%d_rreinit%d_pcg%d%s',...
            cfg.sig_str,cfg.time_jump,cfg.dt,cfg.nt,cfg.first_time,cfg.last_time,cfg.beta,cfg.alpha,cfg.smooth,cfg.dTri,cfg.reinitR,cfg.niter_pcg,scale_str);
        cfg.out_dir=sprintf('%s/%s/%s',cfg.out_dir_pre,cfg.data_tag,cfg.version);
        %%
    case 'gauss' % for test
        cfg.data_tag = 'gauss';
        cfg.data_dir = '';
        cfg.data_name = '';
        cfg.extension = '';
        
        %cfg.anato = '';
        %cfg.sp_mask_opts;
        
        cfg.dataset_name = 'Gauss';
        
        cfg.mask_number = 0;%1;
        switch cfg.mask_number
            case 0 % no mask
                cfg.x_range = 1:20;
                cfg.y_range = 1:20;
                cfg.z_range = 1:20;
            case 1
                cfg.ROI_msk_path = '';
                cfg.x_range = 1:20;
                cfg.y_range = 1:20;
                cfg.z_range = 1:20;
        end
        
        cfg.domain_size = [20,20,20];
        cfg.n1=length(cfg.x_range);
        cfg.n2=length(cfg.y_range);
        cfg.n3=length(cfg.z_range);
        
        cfg.do_resize = 0;
        cfg.size_factor = 1;
        
        cfg.true_size=round(cfg.size_factor*[length(cfg.x_range),length(cfg.y_range),length(cfg.z_range)]);
        
        cfg.reinitR = 0;%0 if do consecutively and 1 if reinitialize rho
        cfg.smooth = 0;%1;
        
        cfg.first_time = 1;
        cfg.time_jump = 1;
        cfg.last_time = 1;
        
        cfg.sigma = 2e-3;
        cfg.sig_str = '2e3';
        cfg.dt = 0.4;
        cfg.nt = 10;
        cfg.alpha = 1000;%100;%1000; % weight for source
        cfg.beta  = 5000;%500;%5000; % weight for matching final image
        cfg.niter_pcg = 60;
        
        cfg.dTri = 1;%3;%1 := 'closed', 3:= 'open'

        cfg.version = sprintf('diff_%s_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_beta_%d_alpha_%d_smooth%d_dtri%d_rreinit%d_pcg%d',...
            cfg.sig_str,cfg.time_jump,cfg.dt,cfg.nt,cfg.first_time,cfg.last_time,cfg.beta,cfg.alpha,cfg.smooth,cfg.dTri,cfg.reinitR,cfg.niter_pcg);
        cfg.out_dir=sprintf('%s/%s/%s',cfg.out_dir_pre,cfg.data_tag,cfg.version);
        %%
    otherwise
        disp('data-name unrecognized')
end
