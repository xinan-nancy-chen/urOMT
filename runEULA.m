function [cfg, Eul] = runEULA(cfg)
% This function runs the Eulerian postprocessing.
addpath('./analyzeFlows','./Inverse','./Sensitivities',genpath('./utilities'));

%%

fprintf('\n =============== Eulerian Post-processing Starts ===============\n')
fprintf('_________________________________________________________\n\n')

tic

date_str = datestr(now,'mmddyy');
paper_fig_str = sprintf('set0%02d',1);


%% get EULA parameters
eulcfg = paramInitEULApar(cfg);

%% set usefule variables:
tag = cfg.data_tag;
cfg.tag = cfg.data_tag;
nt = cfg.nt;
n = cfg.true_size;
ti = cfg.first_time;
tf = cfg.last_time;
tj = cfg.time_jump;
sig_str = cfg.sig_str;

cfg.n = n';
h1 = cfg.spacing.x; h2 = cfg.spacing.y; h3 = cfg.spacing.z;
cfg.h1 = h1.*ones(n(1),1);
cfg.h2 = h2.*ones(n(2),1);
cfg.h3 = h3.*ones(n(3),1);
cfg.dim = length(cfg.n);

switch cfg.dTri
    case 1
        cfg.bc = 'closed';
    case 3
        cfg.bc = 'open';
end
[cfg.Xc,cfg.Yc,cfg.Zc] = getCellCenteredGrid(cfg.h1,cfg.h2,cfg.h3);
cfg.Grad = getCellCenteredGradMatrix({'ccn' 'ccn' 'ccn'},cfg.h1,cfg.h2,cfg.h3);

msk = cfg.msk;

%% initiate
% directory of saved output
outdir = sprintf('EULA_%s_%s/Speed',paper_fig_str,date_str);
cfg.outdir_Eul = outdir;
outdir_long = sprintf('EULA_%s_type_%s_smoothv%s_smoothr%s_%s_%s',...
        cfg.tag,'speed',eulcfg.smoothvSTR,eulcfg.smoothrSTR,paper_fig_str,date_str);

if ~exist(sprintf('%s/%s',cfg.out_dir,outdir),'dir')
    mkdir(sprintf('%s/%s',cfg.out_dir,outdir))
end

fid = fopen(sprintf('%s/%s/%s_record_%s_%s.txt',cfg.out_dir,outdir,cfg.tag,paper_fig_str,date_str),'a+');
fprintf(fid,'%s/%s/%s_record_%s_%s\noutdir-long: %s',cfg.out_dir,outdir,cfg.tag,paper_fig_str,date_str,outdir_long);
title_str = sprintf('============= initiating...\n\nEulerian Analysis %s\n smoothv%s_smoothr%s_%s_%s\n\n',...
    cfg.tag,eulcfg.smoothvSTR,eulcfg.smoothrSTR,paper_fig_str,date_str);
fprintf(title_str)

%
npoints = prod(n);
Ev = NaN(npoints,3,length(ti:tj:tf)*nt);
Es = NaN(npoints,length(ti:tj:tf)*nt);
Eds = NaN(npoints,length(ti:tj:tf)*nt);
Eintp = NaN(npoints,length(ti:tj:tf)*nt+1);
Er = NaN(npoints,length(ti:tj:tf)*nt);

Eintp(:,1) = cfg.vol(1).data(:);

%% running loops
Uall = zeros(3*prod(n),nt*length(ti:tj:tf));
for t1 = ti:tj:tf
    U = load(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,t1,t1+tj,(t1-ti)/tj + 1), 'u');
    Uall(:,(t1-ti)/tj*nt+1:(t1-ti)/tj*nt+nt) = reshape(U.u,[],nt);
end
if eulcfg.smoothv && eulcfg.Svt>0 % smooth velocity field in time space
    fprintf('\n\n============ Smoothing velocity field in time space\n\n')
    Uall = smoothdata(Uall,2,eulcfg.Svt_type,eulcfg.Svt);
end

Rall = zeros(prod(n),nt*length(ti:tj:tf));
for t1 = ti:tj:tf
    R = load(sprintf('%s/r0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,t1,t1+tj,(t1-ti)/tj + 1), 'r');
    Rall(:,(t1-ti)/tj*nt+1:(t1-ti)/tj*nt+nt) = reshape(R.r,[],nt);
end
if eulcfg.smoothr && eulcfg.Srt>0 % smooth relative source in time space
    fprintf('\n\n============ Smoothing relative source in time space\n\n')
    Rall = smoothdata(Rall,2,eulcfg.Srt_type,eulcfg.Srt);
end

tind = 1;
for t1 = ti:tj:tf
    %U = load(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,t1,t1+tj,(t1-ti)/tj + 1), 'u');
    %U = reshape(U.u,[],nt);
    U = Uall(:,(t1-ti)/tj*nt+1:(t1-ti)/tj*nt+nt);
    R = Rall(:,(t1-ti)/tj*nt+1:(t1-ti)/tj*nt+nt);

    if eulcfg.smoothv==1 && eulcfg.Svs>0 % smooth velocity field in spatial space
    for t2 = 1:nt
        u = reshape(U(:,t2),[],3);
        zz = smoothn({reshape(u(:,1),n),reshape(u(:,2),n),reshape(u(:,3),n)},eulcfg.Svs);
        v1 = zz{1};
        v2 = zz{2};
        v3 = zz{3};
        U(:,t2) = [v1(:);v2(:);v3(:)];
    end
    end

    %R = load(sprintf('%s/r0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,t1,t1+tj,(t1-ti)/tj + 1), 'r');
    %R = reshape(R.r,[],nt);
    
    if t1 == ti
        RHO = load(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.tag,ti));
    else
        RHO = load(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,t1-tj,t1,(t1-ti)./tj));
    end
    RHO = RHO.rho_n;
    if isfield(cfg,'chi')
        par = paramInitFunc(cfg.true_size',cfg.nt,cfg.dt,cfg.sigma,cfg.alpha,cfg.beta,cfg.niter_pcg,cfg.dTri,cfg.spacing,cfg.chi((t1-ti)/tj+1).data);
    else
        par = paramInitFunc(cfg.true_size',cfg.nt,cfg.dt,cfg.sigma,cfg.alpha,cfg.beta,cfg.niter_pcg,cfg.dTri,cfg.spacing);
    end
    RHO_t = [RHO, SourceAdvecDiff(RHO,U(:),R(:),nt,par.dt,par)];
    
    Eintp(:,(tind-1)*nt+1+1:tind*nt+1) = RHO_t(:,2:end);%RHO_t(:,1:nt);
    
    for t2 = 1:nt
        TIND = ((t1 - ti)/tj)*nt + t2;
        T = t1+(t2-1)*(tj/nt);
        fprintf('t = %d (t1 = %d, t2 = %d -> T = %.3f)\n',TIND,t1,t2,T);
        
        u = reshape(U(:,t2),[],3); 
        d = reshape(RHO_t(:,t2),n);
        [w2,w1,w3] = gradient(log(d+2*eps)); w1 = w1/cfg.spacing.x; w2 = w2/cfg.spacing.y; w3 = w3/cfg.spacing.z;
        du = cfg.sigma*[w1(:),w2(:),w3(:)];
        Es(:,(tind-1)*nt+t2) = sqrt(u(:,1).^2+u(:,2).^2+u(:,3).^2);
        Eds(:,(tind-1)*nt+t2) = sqrt(du(:,1).^2+du(:,2).^2+du(:,3).^2);
        Er(:,(tind-1)*nt+t2) = R(:,t2);
        Ev(:,1,(tind-1)*nt+t2) = u(:,1);
        Ev(:,2,(tind-1)*nt+t2) = u(:,2);
        Ev(:,3,(tind-1)*nt+t2) = u(:,3);
    end
    tind = tind + 1;
end
Epe = Es./(Eds+2*eps);

T = toc;
fprintf('\n Post Processing --- Elapsed Time: %s\n',datestr(seconds(T),'HH:MM:SS'))

%%
Eul.Es = Es;
Eul.Epe = Epe;
Eul.Er = Er;
Eul.Eintp = Eintp;

end
