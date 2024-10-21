function [cfg,flag] = runUROMT(cfg)
% This function runs the urOMT algorithm.


%% %%%%%%%%%%%%%%
if ~exist(cfg.out_dir,'dir')
    mkdir(cfg.out_dir)
end

reInitializeU = 1;

fname=sprintf('%s/record.txt',cfg.out_dir);
if ~exist(sprintf('%s/record.txt',cfg.out_dir),'file')
    csvwrite_with_headers(fname,[0 0 0 0 0 0 0 0 0 0 0 0],{'time-ind','ti','tf','gamma','gamma1','gamma2','gamma3','max(u)','min(r)','max(r)','toc','alpha'});
end

rho_n = cfg.vol(1).data(:);

if ~exist(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.data_tag,cfg.first_time),'file')
    save_rhon(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.data_tag,cfg.first_time),rho_n);
end
fprintf('\n =============== urOMT Starts ===============\n')
fprintf('______________________________________________\n\n')
fprintf(' tag:\t\t%s\n sigma:\t\t%.4f\n alpha:\t\t%.4f\n beta:\t\t%.4f\n nt:\t\t%d\n dt:\t\t%.2f\n pcg:\t\t%d\n',cfg.data_tag,cfg.sigma,cfg.alpha,cfg.beta,cfg.nt,cfg.dt,cfg.niter_pcg)
fprintf(' size:\t\t[%d,%d,%d]\n do_resize:\t%d\n resize_factor:\t%.2f\n start frame:\t%d\n end frame:\t%d\n frame jump:\t%d\n if use chi:\t%d\n\n\n',cfg.true_size(1),cfg.true_size(2),cfg.true_size(3),cfg.do_resize,cfg.size_factor,cfg.first_time,cfg.last_time+cfg.time_jump,cfg.time_jump,isfield(cfg,'chi'))

% main loop
t1 = tic;
for tind = 1:length(cfg.first_time:cfg.time_jump:cfg.last_time)
    fprintf('tind = %d\n',tind)
    t2 = tic;
    
    if exist(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.data_tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),'file')==2 && exist(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.data_tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),'file') == 2
        rho_n = load_rhon(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.data_tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind));
        u = load_un(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.data_tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind));
        u = reshape(u,[],cfg.nt);
        continue
    end
    
    %source density
    if cfg.reinitR
        rho_0 = cfg.vol(tind).data;
        rho_0 = rho_0(:);
    else
        rho_0 = rho_n(:);
    end

    %target density
    rho_N = cfg.vol(tind+1).data;
    
    % set parameters
    if isfield(cfg,'chi')
        par = paramInitFunc(cfg.true_size',cfg.nt,cfg.dt,cfg.sigma,cfg.alpha,cfg.beta,cfg.niter_pcg,cfg.dTri,cfg.spacing,cfg.chi(tind).data);
    else
        par = paramInitFunc(cfg.true_size',cfg.nt,cfg.dt,cfg.sigma,cfg.alpha,cfg.beta,cfg.niter_pcg,cfg.dTri,cfg.spacing); 
    end
    
    par.drho0     = rho_0(:);
    par.drhoN     = rho_N(:);
    
    % initial guess for u &r:
    if tind == 1 || reInitializeU
        u = zeros(par.dim*prod(par.n),par.nt);
        r = zeros(prod(par.n),par.nt);
    end
    
    %% Descent for u & r
    fprintf('\n ============================== Descent on u & r ==============================\n')
    fprintf('__________________________________________________________________________________________________\n\n')
    fprintf('i.lsiter\t   Gamma(ttl)  \tGamma1(kint)    \tGamma2(src)    \tGamma3(fit)    \t descent output\n')
    fprintf('__________        ____________     ____________     ____________     ____________      ___________________\n')
    [u,r,GAMMA,~] = GNblock_ur(u,r,par.nt,par.dt,par,cfg.data_tag);
    
    Rho = SourceAdvecDiff(par.drho0,u,r,par.nt,par.dt,par);
    rho_n = Rho(:,end);
    btoc = toc(t2);
    
    dlmwrite(fname,[tind,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,GAMMA.Gamma,GAMMA.Gamma1,GAMMA.Gamma2,GAMMA.Gamma3,max(u(:)),min(r(:)),max(r(:)),btoc,par.alpha],'-append');
    
    save_un(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.data_tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),u);
    save_rn(sprintf('%s/r0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.data_tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),r);
    save_rhon(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.data_tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),rho_n);
    
    fprintf('tind = %d, max(u) = %5.4f, min(r) = %5.4f, max(r) = %5.4f\n',tind,max(u(:)),min(r(:)),max(r(:)));
end
Btoc = toc(t1);
dlmwrite(fname,[0,0,0,0,0,0,0,0,0,0,Btoc],'-append');
fprintf('\n =============== urOMT Ends ===============\n')
fprintf('\n Elapsed Time: %s\n',datestr(seconds(Btoc),'DD:HH:MM:SS'))

flag = 1;

end