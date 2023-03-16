function [cfg, Lag] = runGLAD(cfg)
% This function runs the Lagrangian postprocessing.
addpath('./analyzeFlows','./Inverse','./Sensitivities',genpath('./utilities'));

%%

fprintf('\n =============== Lagrangian Post-processing Starts ===============\n')
fprintf('_________________________________________________________\n\n')

tic

date_str = datestr(now,'mmddyy');
paper_fig_str = sprintf('set0%02d',1);


%% get GLAD parameters
glacfg = paramInitGLADpar(cfg);

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

%%
if glacfg.do_sp
    data_max = zeros(size(cfg.vol(1).data));
    for l = 1:length(cfg.vol)
        data_max = max(data_max, cfg.vol(l).data);
    end
    mind = find((msk>0) & (data_max>glacfg.sp_thresh));
    glacfg.do_sp_str = sprintf('_data_min_%d',glacfg.sp_thresh);
else
    mind = find(msk>0);
    glacfg.do_sp_str = '';
end

[syind,sxind,szind] = ind2sub(n,mind); %find indices of all voxels inside the ROI-SP
switch glacfg.spType
    case 'ordered'
        sy = syind(1:glacfg.spfs:end);
        sx = sxind(1:glacfg.spfs:end);
        sz = szind(1:glacfg.spfs:end);
    case 'uniform'
        mskvol = sum(msk(:)); %volume of mask used to select start points
        NSP = round(glacfg.spPerc*mskvol/100);
        if ~exist(sprintf('%s_%s%d_spPerc%d_nsp%d.mat',tag,glacfg.spMsk_name,glacfg.spMsk_ind,glacfg.spPerc,NSP),'file')
            [spIND,spINDid] = datasample(mind,NSP,'Replace',false);
            save(sprintf('%s_%s%d_spPerc%d_nsp%d.mat',tag,glacfg.spMsk_name,glacfg.spMsk_ind,glacfg.spPerc,NSP),'spIND');
        else
            load(sprintf('%s_%s%d_spPerc%d_nsp%d.mat',tag,glacfg.spMsk_name,glacfg.spMsk_ind,glacfg.spPerc,NSP));
        end
        [sy,sx,sz] = ind2sub(n,sort(spIND,'ascend'));
end



%initialize streamlines:
nsp = length(sx);
%convert from matlab grid to cell-centered grid:
s1 = (sy-0.5).*h1; %i/y-axis
s2 = (sx-0.5).*h2; %j/x-axis
s3 = (sz-0.5).*h3; %k/z-axis
sp_123 = [s1,s2,s3];


pcur = sp_123; %current point i.e. list of current location in each streamline that hasn't been terminated
npoints = length(pcur); %keep track of the # of streamlines that have not yet been terminated 

SL = cell(1,nsp);
SPD_SL = cell(1,nsp);
DSPD_SL = cell(1,nsp);
AUGSPD_SL = cell(1,nsp);
RHO_SL = cell(1,nsp);
PE_SL = cell(1,nsp);
FLX_SL = cell(1,nsp);
AUGFLX_SL = cell(1,nsp);
R_SL = cell(1,nsp);

%initiate pathlines
pl = NaN(npoints,3,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
plspd = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
pldspd = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
plaugspd = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
plrho = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
plflx = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
plaugflx = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
plr = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);

step = 1;
xt = pcur; %current position in rcl orientation
xt = max([h1*0.5,h2*0.5,h3*0.5],min(xt,[h1*(n(1)-.5001),h2*(n(2)-.5001),h3*(n(3)-.5001)]));

%
Uall = zeros(3*prod(n),nt*length(ti:tj:tf));
for t1 = ti:tj:tf
    U = load(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,t1,t1+tj,(t1-ti)/tj + 1), 'u');
    Uall(:,(t1-ti)/tj*nt+1:(t1-ti)/tj*nt+nt) = reshape(U.u,[],nt);
end
if glacfg.smoothv && glacfg.Svt>0 % smooth velocity field in time space
    Uall = cell2mat(smoothn(num2cell(Uall',1),glacfg.Svt));
    Uall = Uall';
end

% running pathlines

for t1 = ti:tj:tf
    U = Uall(:,(t1-ti)/tj*nt+1:(t1-ti)/tj*nt+nt);
    R = load(sprintf('%s/r0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,t1,t1+tj,(t1-ti)/tj + 1), 'r');
    R = reshape(R.r,[],nt);
    
    if glacfg.smoothv==1 && glacfg.Svs>0 % smooth velocity field in spatial space
        for t2 = 1:nt
            u = reshape(U(:,t2),[],3);
            zz = smoothn({reshape(u(:,1),n),reshape(u(:,2),n),reshape(u(:,3),n)},glacfg.Svs);
            v1 = zz{1};
            v2 = zz{2};
            v3 = zz{3};
            U(:,t2) = [v1(:);v2(:);v3(:)];
        end
    end
    if strcmp(glacfg.RD,'R')
        if t1 == ti
            RHO = load(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.tag,ti));
        else
            RHO = load(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,t1-tj,t1,(t1-ti)./tj));
        end
        RHO = RHO.rho_n;
        if nt > 1
            if isfield(cfg,'chi')
                par = paramInitFunc(cfg.true_size',cfg.nt,cfg.dt,cfg.sigma,cfg.alpha,cfg.beta,cfg.niter_pcg,cfg.dTri,cfg.spacing,cfg.chi((t1-ti)/tj+1).data);
            else
                par = paramInitFunc(cfg.true_size',cfg.nt,cfg.dt,cfg.sigma,cfg.alpha,cfg.beta,cfg.niter_pcg,cfg.dTri,cfg.spacing);
            end
            RHO_t = [RHO, SourceAdvecDiff(RHO,U(:),R(:),nt,par.dt,par)];
        else
            RHO_t = RHO;
        end
    end
    
    for t2 = 1:nt
        TIND = ((t1 - ti)/tj)*nt + t2;
        T = t1+(t2-1)*(tj/nt);
        fprintf('t = %d (t1 = %d, t2 = %d -> T = %.3f)\n',TIND,t1,t2,T);
        
        switch glacfg.RD
            case 'D'
                d = cfg.vol(round((t1 - ti)/tj + 1 + (t2-1)/nt)).data;
            case 'R'
                d = reshape(RHO_t(:,t2),n);
        end
        
        if glacfg.minIm0
            d = d-min(d(:));
        end
        
        %compute streamlines
        u = reshape(U(:,t2),[],3);
        v1 = reshape(u(:,1),n);
        v2 = reshape(u(:,2),n);
        v3 = reshape(u(:,3),n);
        
        Rr = reshape(R(:,t2),n);
        
        %add eps to d so can take log(d) and not log(0)
        [w2,w1,w3] = gradient(log(d+2*eps)); w1 = w1/cfg.spacing.x; w2 = w2/cfg.spacing.y; w3 = w3/cfg.spacing.z;
        u1 = v1 - cfg.sigma.*w1;
        u2 = v2 - cfg.sigma.*w2;
        u3 = v3 - cfg.sigma.*w3;
        
        du = cfg.sigma*[w1(:),w2(:),w3(:)];
        
        switch glacfg.flw_type
            case 'vel'
                a1=u1;
                a2=u2;
                a3=u3;
            case 'flw'
                a1=u1.*d;
                a2=u2.*d;
                a3=u3.*d;
        end

        speed = reshape(sqrt(sum(u.^2,2)),n);
        dspeed = reshape(sqrt(sum(du.^2,2)),n);
        speedAug = sqrt(u1.^2 + u2.^2 + u3.^2);
        img_flow = speed.*d;
        img_flowAug = speedAug.*d;
        
        switch glacfg.cutoff_str
            case 'min'
                conf.conc = 1;
                conf.speed = 1;
            case 'max'
                conf.conc = mean(d(d>0)) + std(d(d>0));
                conf.speed = mean(speed(speed>0)) + std(speed(speed>0));              
            case 'mean'
                conf.conc = mean(d(d>0));
                conf.speed = mean(speed(speed>0));
        end
        
        % update first step of density and speed
        if step == 1
            pl(:,:,1) = xt;
            plspd(:,1) = speed(sub2ind(n,sy,sx,sz));
            pldspd(:,1) = dspeed(sub2ind(n,sy,sx,sz));
            plaugspd(:,1) = speedAug(sub2ind(n,sy,sx,sz));
            plrho(:,1) = d(sub2ind(n,sy,sx,sz));
            plflx(:,1) = img_flow(sub2ind(n,sy,sx,sz));
            plaugflx(:,1) = img_flowAug(sub2ind(n,sy,sx,sz));
            plr(:,1) = Rr(sub2ind(n,sy,sx,sz));
        end
        
        %vector field to be integrated in order to compute streamlines
        V = [a1(:),a2(:),a3(:)];
        V(msk==0,:) = 0;
        u(msk==0,:) = 0;
        du(msk==0,:) = 0;
        
        for Estep = 1:glacfg.nEulStep
            step = step + 1;
            V_interp = interp_vel(V,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            u_interp = interp_vel(u,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            du_interp = interp_vel(du,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            
            conc_interp = interpF(d,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            r_interp = interpF(Rr,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            spd_interp = sqrt(sum(u_interp.^2,2));
            dspd_interp = sqrt(sum(du_interp.^2,2));
            spdAug_interp = sqrt(sum(V_interp.^2,2));
            flx_interp = interpF(img_flow,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            flxaug_interp = interpF(img_flowAug,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            
            conc_conf_interp = conc_interp./conf.conc;
            spd_conf_interp = spd_interp./conf.speed;
            
            thrsh_ind = find(conc_conf_interp > glacfg.thresholds.conc & spd_conf_interp > glacfg.thresholds.speed);

            if isempty(thrsh_ind)
                break
            else
                switch glacfg.XT
                    case'T'
                        xt(thrsh_ind,:) = xt(thrsh_ind,:) + glacfg.mdt.*V_interp(thrsh_ind,:);%so keep all current locations in case can take a step later when conc gets there
                    case 'X'
                        xt(thrsh_ind,:) = xt(thrsh_ind,:) + glacfg.mdt.*(V_interp(thrsh_ind,:)./spd_interp);%so keep all current locations in case can take a step later when conc gets there
                end
                %make sure it stays in bounds:
                xt = max([h1*0.5,h2*0.5,h3*0.5],min(xt,[h1*(n(1)-.5001),h2*(n(2)-.5001),h3*(n(3)-.5001)]));
                pl(thrsh_ind,:,step) = xt(thrsh_ind,:);
                plspd(thrsh_ind,step) = spd_interp(thrsh_ind);
                pldspd(thrsh_ind,step) = dspd_interp(thrsh_ind);
                plaugspd(thrsh_ind,step) = spdAug_interp(thrsh_ind);
                plrho(thrsh_ind,step)  = conc_interp(thrsh_ind);
                plflx(thrsh_ind,step) = flx_interp(thrsh_ind);
                plaugflx(thrsh_ind,step) = flxaug_interp(thrsh_ind);
                plr(thrsh_ind,step)  = r_interp(thrsh_ind);
                
            end
        end
    end
end
clear Uall

% select qualified pathlines
pli = 0; %counter for added pathlines
for sli = 1:npoints
    pl_cur = squeeze(pl(sli,:,:))';
    aind = any(~isnan(pl_cur),2); %1 if row is not NaN, 0 if row is NaN
    pl_cur = pl_cur(aind,:);
    %check that unique sline has more than 1 point (remove not-move coordinates)
    [pl_cur,ia,~] = unique(pl_cur,'rows','stable');
    if size(pl_cur,1)>glacfg.pln
        pli = pli + 1;
        if glacfg.smoothp == 1 && size(pl_cur,1)>3 % smooth pathline in spatial space
            zz = smoothn({pl_cur(2:end-1,1),pl_cur(2:end-1,2),pl_cur(2:end-1,3)},glacfg.smpTol);
            pl_cur = [[pl_cur(1,1);zz{1};pl_cur(end,1)],[pl_cur(1,2);zz{2};pl_cur(end,2)],[pl_cur(1,3);zz{3};pl_cur(end,3)]];
        end
        SL{pli} = pl_cur;
        
        pls_cur = plspd(sli,aind)';
        plds_cur = pldspd(sli,aind)';
        plr_cur = plrho(sli,aind)';
        plas_cur = plaugspd(sli,aind)';
        plflx_cur = plflx(sli,aind)';
        plaugflx_cur = plaugflx(sli,aind)';
        plRr_cur = plr(sli,aind)';
        
        SPD_SL{pli} = pls_cur(ia);
        DSPD_SL{pli} = plds_cur(ia);
        AUGSPD_SL{pli} = plas_cur(ia);
        RHO_SL{pli} = plr_cur(ia);
        PE_SL{pli} = pls_cur(ia)./(plds_cur(ia)+2*eps);
        FLX_SL{pli} = plflx_cur(ia);
        AUGFLX_SL{pli} = plaugflx_cur(ia);
        R_SL{pli} = plRr_cur(ia);
    end
end

%remove empty cell spaces from pathlines that were removed:
SL = SL(1:pli);
SPD_SL = SPD_SL(1:pli);
DSPD_SL = DSPD_SL(1:pli);
AUGSPD_SL = AUGSPD_SL(1:pli);
RHO_SL = RHO_SL(1:pli);
PE_SL = PE_SL(1:pli);
FLX_SL = FLX_SL(1:pli);
AUGFLX_SL = AUGFLX_SL(1:pli);
R_SL = R_SL(1:pli);

sl_euc = cellfun(@(x) sqrt(sum((x(end,:)-x(1,:)).^2)),SL);% getcell array with euclidean length of sl between first and last point
%figure, subplot(2,2,1),histogram(sl_euc),title('sl-euc'),axis tight; subplot(2,2,2),histogram(cat(1, R_SL{:})),title('R-relative source'),axis tight; subplot(2,2,3),histogram(cat(1, SPD_SL{:})),title('S-spd'),axis tight; subplot(2,2,4),histogram(cat(1, PE_SL{:})),title('Pe-Peclet'),axis tight, xlim([0,1000]); set(gcf,'Position',[560   375   747   573]);

SL2 = SL(sl_euc>glacfg.sl_tol);%only keep streamlines whose Euclidean length between first and last points is larger than the threshold
sstream2 = SPD_SL(sl_euc>glacfg.sl_tol);
dsstream2 = DSPD_SL(sl_euc>glacfg.sl_tol);
asstream2 = AUGSPD_SL(sl_euc>glacfg.sl_tol);
rstream2 = RHO_SL(sl_euc>glacfg.sl_tol);
pestream2 = PE_SL(sl_euc>glacfg.sl_tol);
fstream2 = FLX_SL(sl_euc>glacfg.sl_tol);
afstream2 = AUGFLX_SL(sl_euc>glacfg.sl_tol);
Rstream2 = R_SL(sl_euc>glacfg.sl_tol);

fprintf(' # of start points = %d\n # of effective pathlines after pathline-number (pln) threshold = %d \n # of effective pathlines after Euclidean dist (sl_tol) threshold = %d\n',npoints,pli,length(SL2))
pl_cur = cellfun(@(x) x(:,[2,1,3]),SL2,'UniformOutput',false);

% s

outdir = sprintf('LPPA_%s_%s/Speed',paper_fig_str,date_str);
cfg.outdir_Lag_s = outdir;
outdir_long = sprintf('LPPA_%s_type_%s_%s_smoothv%s_smoothp%s_img%s_mdt%d%s_%s_nEul%d_cutoffStr_%s_concThresh%5.4f_spdThresh%5.4f_minIm0_%d_%s_slTol%d_sfs%d_%s_diff%s_tj%d_nt%d%s_%s_%s',...
        cfg.tag,'speedmap',glacfg.flw_type,glacfg.smoothvSTR,glacfg.smoothpSTR,glacfg.RD,glacfg.mdt,glacfg.XT,glacfg.spMSK_str,glacfg.nEulStep,glacfg.cutoff_str,glacfg.thresholds.conc,glacfg.thresholds.speed,glacfg.minIm0,glacfg.spSTR,glacfg.sl_tol,glacfg.sfs,glacfg.QB_str,sig_str,tj,nt,glacfg.do_sp_str,paper_fig_str,date_str);

if ~exist(sprintf('%s/%s',cfg.out_dir,outdir),'dir')
    mkdir(sprintf('%s/%s',cfg.out_dir,outdir))
end

%saveas(gcf,sprintf('%s/%s/%s_Lagresults_stat_E%02d_%02d_%s_%s.png',cfg.out_dir,outdir,cfg.tag,ti,tf+tj,paper_fig_str,date_str));

fid = fopen(sprintf('%s/%s/%s_record_%s_%s.txt',cfg.out_dir,outdir,cfg.tag,paper_fig_str,date_str),'a+');
fprintf(fid,'%s/%s/%s_record_%s_%s\noutdir-long: %s',cfg.out_dir,outdir,cfg.tag,paper_fig_str,date_str,outdir_long);
title_str = sprintf('============= initiating...\n\nLagrangian-Pathline (%s data, affSmooth = %d), \nanalysis type = %s\n \nflw_type = %s, smoothv = %s, smoothp = %s, img = %s, mdt = %d(%s), %s, nEulStep= %d, \ncutoffStr = %s, concThresh = %5.4f, spdThresh = %5.4f, minIm0 = %d, %s, slTol = %d, sfs = %d\n %s, \ndiff = %s, tj = %d, nt= %d%s_%s_%s\n\n',...
    cfg.tag,cfg.smooth,'speedmap',glacfg.flw_type,glacfg.smoothvSTR,glacfg.smoothpSTR,glacfg.RD,glacfg.mdt,glacfg.XT,glacfg.spMSK_str,glacfg.nEulStep,glacfg.cutoff_str,glacfg.thresholds.conc,glacfg.thresholds.speed,glacfg.minIm0,glacfg.spSTRtitle,glacfg.sl_tol,glacfg.sfs,glacfg.QB_str,sig_str,tj,nt,glacfg.do_sp_str,paper_fig_str,date_str);
fprintf(title_str)

fprintf('getting speed map...\n')
%
% initialize temporary masks:
s = zeros(n);%speed
ds = zeros(n);%diffusive speed
r = zeros(n);%relative source
s_full = zeros(n);%full smooth speed
ds_full = zeros(n);%full smooth diffusive speed
%as_full = zeros(n);%full smooth augmented speed
r_full = zeros(n);%full relative source
scount = zeros(n);%full count
%f_full = zeros(n);%full smooth flux
%af_full = zeros(n);%full smooth augmented flux

if exist('pl_cur.mat','file')
    delete('pl_cur.mat');
end
if exist('pli_array.mat','file')
    delete('pli_array.mat');
end

if glacfg.do_QB
    % average map
    for k = 1:length(pl_cur)
        % streamline cluster mask:
        slines_tmp = pl_cur(k);
        spdlines_tmp = sstream2(k);
        dspdlines_tmp = dsstream2(k);

        n_slines = size(slines_tmp,2);
        for ind_line = 1:n_slines
            %convert back to MATLAB grid
            sline = round((slines_tmp{ind_line}./[h1,h2,h3])+0.5);
            [sline,ia,~] = unique(sline, 'rows', 'stable');
            ssl = spdlines_tmp{ind_line}(ia);
            dssl = dspdlines_tmp{ind_line}(ia);

            subs_1 = sline(:,2);
            subs_2 = sline(:,1);
            subs_3 = sline(:,3);
            subs_1(subs_1 < 1) = 1;
            subs_2(subs_2 < 1) = 1;
            subs_3(subs_3 < 1) = 1;
            subs_1(subs_1 > n(1)) = n(1);
            subs_2(subs_2 > n(2)) = n(2);
            subs_3(subs_3 > n(3)) = n(3);
            inds = sub2ind(n, subs_1, subs_2, subs_3);

            s_full(inds) = ssl + s_full(inds);
            ds_full(inds) = dssl + ds_full(inds);
            scount(inds) = scount(inds) + 1;
        end
    end
    Pe_full = s_full./ds_full;%s_full./(ds_full+eps);
    s_full(s_full>0) = s_full(s_full>0)./scount(scount>0);
    
    % non-average map
    pl_cur = pl_cur(1:glacfg.sfs:end);
    sstream3 = sstream2(1:glacfg.sfs:end);
    dsstream3 = dsstream2(1:glacfg.sfs:end);
    save('pl_cur.mat','pl_cur');

    fprintf('%s: waiting to run QB.py\n',cfg.tag)
    pause
    fprintf('continue\n')

    load(fullfile(cd,'pli_array.mat'));
    % only want clusters with more than tol # of sls:
    clus_length = cellfun('size',pli_array,2)'; fprintf(' # of original clusters = %d\n',length(clus_length))
    sli = pli_array(clus_length>glacfg.clus_tol);
    nclus = size(sli,2);%# of clusters that made the cutoff
    clus_length = cellfun('size',sli,2)'; fprintf(' # of clusters after cluster-length (clus_tol) threshold = %d\n',length(clus_length))
    [~,clus_size_rank] = sort(clus_length,'descend');
    %only keep largest N clusters where N = clus_cutoff:
    if glacfg.clus_cutoff>0
        if nclus > glacfg.clus_cutoff
            ind_tmp = zeros(nclus,1);
            ind_tmp(clus_size_rank(1:glacfg.clus_cutoff)) = 1;
            sli = sli(ind_tmp==1);

            nclus = size(sli,2);

            clus_length = cellfun('size',sli,2)';
            [~,clus_size_rank] = sort(clus_length);
        end
    end
    fprintf(' # of clusters after max-cluster-number (clus_cutoff) threshold = %d\n',nclus)

    %getting clustered pathline start points
    for k = 1:nclus
        % streamline cluster mask:
        slines_tmp = pl_cur([sli{k}]+1);
        spdlines_tmp = sstream3([sli{k}]+1);
        dspdlines_tmp = dsstream3([sli{k}]+1);

        n_slines = size(slines_tmp,2);
        for ind_line = 1:n_slines
            %convert back to MATLAB grid
            sline = round((slines_tmp{ind_line}./[h1,h2,h3])+0.5);
            [sline,ia,~] = unique(sline, 'rows', 'stable');
            ssl = spdlines_tmp{ind_line}(ia);
            dssl = dspdlines_tmp{ind_line}(ia);

            subs_1 = sline(:,2);
            subs_2 = sline(:,1);
            subs_3 = sline(:,3);
            subs_1(subs_1 < 1) = 1;
            subs_2(subs_2 < 1) = 1;
            subs_3(subs_3 < 1) = 1;
            subs_1(subs_1 > n(1)) = n(1);
            subs_2(subs_2 > n(2)) = n(2);
            subs_3(subs_3 > n(3)) = n(3);
            inds = sub2ind(n, subs_1, subs_2, subs_3);

            s(inds) = ssl;
            ds(inds) = dssl;
            %s_full(inds) = ssl + s_full(inds);
            %ds_full(inds) = dssl + ds_full(inds);
            %scount(inds) = scount(inds) + 1;
        end
    end
    %Pe_full = s_full./ds_full;
    %s_full(s_full>0) = s_full(s_full>0)./scount(scount>0);
    %ds_full(ds_full>0) = ds_full(ds_full>0)./scount(scount>0);

else % No QuickBundle
    % non-average map
    for k = 1:glacfg.sfs:length(pl_cur) 
        % streamline cluster mask:
        slines_tmp = pl_cur(k);
        spdlines_tmp = sstream2(k);
        dspdlines_tmp = dsstream2(k);
        Rlines_tmp = Rstream2(k);

        n_slines = size(slines_tmp,2);
        for ind_line = 1:n_slines
            %convert back to MATLAB grid
            sline = round((slines_tmp{ind_line}./[h1,h2,h3])+0.5);
            [sline,ia,~] = unique(sline, 'rows', 'stable');
            ssl = spdlines_tmp{ind_line}(ia);
            dssl = dspdlines_tmp{ind_line}(ia);
            Rsl = Rlines_tmp{ind_line}(ia);

            subs_1 = sline(:,2);
            subs_2 = sline(:,1);
            subs_3 = sline(:,3);
            subs_1(subs_1 < 1) = 1;
            subs_2(subs_2 < 1) = 1;
            subs_3(subs_3 < 1) = 1;
            subs_1(subs_1 > n(1)) = n(1);
            subs_2(subs_2 > n(2)) = n(2);
            subs_3(subs_3 > n(3)) = n(3);
            inds = sub2ind(n, subs_1, subs_2, subs_3);

            s(inds) = ssl;
            ds(inds) = dssl;
            r(inds) = Rsl;
        end
    end
    
    % averaged map
    for k = 1:length(pl_cur)
        % streamline cluster mask:
        slines_tmp = pl_cur(k);
        spdlines_tmp = sstream2(k);
        dspdlines_tmp = dsstream2(k);
        augspdlines_tmp = asstream2(k);
        fluxlines_tmp = fstream2(k);
        afluxlines_tmp = afstream2(k);
        Rlines_tmp = Rstream2(k);
        
        n_slines = size(slines_tmp,2);
        for ind_line = 1:n_slines
            %convert back to MATLAB grid
            sline = round((slines_tmp{ind_line}./[h1,h2,h3])+0.5);
            [sline,ia,~] = unique(sline, 'rows', 'stable');
            ssl = spdlines_tmp{ind_line}(ia);
            dssl = dspdlines_tmp{ind_line}(ia);
            augssl = augspdlines_tmp{ind_line}(ia);
            fluxsl = fluxlines_tmp{ind_line}(ia);
            afluxsl = afluxlines_tmp{ind_line}(ia);
            Rsl = Rlines_tmp{ind_line}(ia);
            
            subs_1 = sline(:,2);
            subs_2 = sline(:,1);
            subs_3 = sline(:,3);
            subs_1(subs_1 < 1) = 1;
            subs_2(subs_2 < 1) = 1;
            subs_3(subs_3 < 1) = 1;
            subs_1(subs_1 > n(1)) = n(1);
            subs_2(subs_2 > n(2)) = n(2);
            subs_3(subs_3 > n(3)) = n(3);
            inds = sub2ind(n, subs_1, subs_2, subs_3);

            s_full(inds) = ssl + s_full(inds);
            ds_full(inds) = dssl + ds_full(inds);
            r_full(inds) = Rsl + r_full(inds);
            scount(inds) = scount(inds) + 1;
        end

    end
    Pe_full = s_full./ds_full;%s_full./(ds_full+eps);
    s_full(scount>0) = s_full(scount>0)./scount(scount>0);
    r_full(scount>0) = r_full(scount>0)./scount(scount>0);

end


if glacfg.do_masked
    s(~msk) = 0;
    ds(~msk) = 0;
    r(~msk) = 0;
    s_full(~msk) = 0;
    ds_full(~msk) = 0;
    Pe_full(~msk) = 0;
    r_full(~msk) = 0;
end


Pe = s./ds; %s./(ds+eps);
Pe(~msk) = 0;
fprintf('tag = %s, detect %d Pe == Nan\n',tag,sum(isnan(Pe(:))));
Pe(isnan(Pe)) = 0;
fprintf('tag = %s, detect %d Pe == Inf\n',tag,length(find(Pe==Inf)));
Pe(Pe==Inf) = 0;%max(Pe(Pe~=Inf));%0;

fprintf('tag = %s, detect %d Pe_full == Nan\n',tag,sum(isnan(Pe_full(:))));
Pe_full(isnan(Pe_full)) = 0;
fprintf('tag = %s, detect %d Pe_full == Inf\n',tag,length(find(Pe_full==Inf)));
Pe_full(Pe_full==Inf) = 0;%max(Pe_full(Pe_full~=Inf));%0;

fprintf('length(s_full(s_full>0)) = %d, length(ds_full(ds_full>0)) = %d, length(Pe_full(Pe_full>0)) = %d\n',length(s_full(s_full>0)),length(ds_full(ds_full>0)),length(Pe_full(Pe_full>0)))

% v

outversion = sprintf('%s_%s',paper_fig_str,date_str);

outdir = sprintf('LPPA_%s_%s/Pathlines',paper_fig_str,date_str);
cfg.outdir_Lag_v = outdir;
outdir_long = sprintf('LPPA_%s_type_%s_%s_smoothv%s_smoothp%s_img%s_mdt%d%s_%s_nEul%d_cutoffStr_%s_concThresh%5.4f_spdThresh%5.4f_minIm0_%d_%s_slTol%d_vfs%d_diff%s_tj%d_nt%d%s_%s_%s',...
        cfg.tag,'vectors',glacfg.flw_type,glacfg.smoothvSTR,glacfg.smoothpSTR,glacfg.RD,glacfg.mdt,glacfg.XT,glacfg.spMSK_str,glacfg.nEulStep,glacfg.cutoff_str,glacfg.thresholds.conc,glacfg.thresholds.speed,glacfg.minIm0,glacfg.spSTR,glacfg.sl_tol,glacfg.vfs,sig_str,tj,nt,glacfg.do_sp_str,paper_fig_str,date_str);
if ~exist(sprintf('%s/%s',cfg.out_dir,outdir),'dir')
    mkdir(sprintf('%s/%s',cfg.out_dir,outdir))
end

fid = fopen(sprintf('%s/%s/%s_record_%s_%s.txt',cfg.out_dir,outdir,cfg.tag,paper_fig_str,date_str),'a+');
fprintf(fid,'%s/%s/%s_record_%s_%s\noutdir-long: %s',cfg.out_dir,outdir,cfg.tag,paper_fig_str,date_str,outdir_long);
title_str = sprintf('============= \n\nLagrangian-Pathline (%s data, affSmooth = %d), \nanalysis type = %s\n \nflw_type = %s, smoothv = %s, smoothp = %s, img = %s, mdt = %d(%s), %s, nEulStep= %d, \ncutoffStr = %s, concThresh = %5.4f, spdThresh = %5.4f, minIm0 = %d, %s, slTol = %d, vfs = %d, \ndiff = %s, tj = %d, nt = %d%s_%s_%s\n\n',...
    cfg.tag,cfg.smooth,'vectors',glacfg.flw_type,glacfg.smoothvSTR,glacfg.smoothpSTR,glacfg.RD,glacfg.mdt,glacfg.XT,glacfg.spMSK_str,glacfg.nEulStep,glacfg.cutoff_str,glacfg.thresholds.conc,glacfg.thresholds.speed,glacfg.minIm0,glacfg.spSTRtitle,glacfg.sl_tol,glacfg.vfs,sig_str,tj,nt,glacfg.do_sp_str,paper_fig_str,date_str);
fprintf(title_str)

fprintf('getting pathlines...\n')
fprintf('Total original %d pathlines\n',length(SL2)); 

%convert from cell-centered grid to matlab grid:
SL = cellfun(@(x) [x(:,1)./h1+0.5,x(:,2)./h2+0.5,x(:,3)./h3+0.5],SL2,'UniformOutput',false); 
clear SL2;

dispnor = cellfun(@(x) (x(end,:)-x(1,:))/norm(x(end,:)-x(1,:)),SL,'UniformOutput',false); 
disp = cellfun(@(x) x(end,:)-x(1,:),SL,'UniformOutput',false); 
startp = cellfun(@(x) x(1,:),SL,'UniformOutput',false); 
endp = cellfun(@(x) x(end,:),SL,'UniformOutput',false); 

PATH.NPtsInPath = cellfun(@(x) size(x,1),SL); % number of points in each pathline
PATH.LengthInPath = cellfun(@(x) sum(sqrt(sum(diff(x).^2,2))),SL); % total length of path in each pathline
PATH.disp = reshape([disp{:}]',3,[])'; % displacement field
PATH.dispnor = reshape([dispnor{:}]',3,[])'; % normalized displacement field
PATH.startp = round(reshape([startp{:}]',3,[])'); % startp points
PATH.endp = reshape([endp{:}]',3,[])'; % endp points
PATH.displen = sqrt(PATH.disp(:,1).^2+PATH.disp(:,2).^2+PATH.disp(:,3).^2); % displacement length in each pathline
%%%%%%%%%%%%here!
if isfield(cfg,'anato')
    anato = load_rhon(cfg.anato); 
    anato = anato(cfg.x_range,cfg.y_range,cfg.z_range);

    if cfg.do_resize
       anato = resizeMatrix(anato,round(cfg.size_factor.*size(anato)),'linear');
    end
    anato(~msk) = 0;
end

indb = find(msk(sub2ind(n,PATH.startp(:,1),round(PATH.startp(:,2)),PATH.startp(:,3)))==1); % index of those starting within brain
PATH.ind_msk = indb;

% save to vtk
%
SL_tmp = SL(indb);
PE_tmp = pestream2(indb);
S_tmp = sstream2(indb);
R_tmp = Rstream2(indb);
vtkwrite_pathlines(sprintf('%s/%s/%s_pathlines_sltol_%.1f_jp_%d_%s.vtk',cfg.out_dir,outdir,tag,glacfg.sl_tol,glacfg.jp,outversion),'polydata','lines',SL_tmp(1:glacfg.jp:end));
vtkwrite_spdlines(sprintf('%s/%s/%s_Pelines_sltol_%.1f_jp_%d_%s.vtk',cfg.out_dir,outdir,tag,glacfg.sl_tol,glacfg.jp,outversion),'polydata','lines',SL_tmp(1:glacfg.jp:end),PE_tmp(1:glacfg.jp:end));
vtkwrite_spdlines(sprintf('%s/%s/%s_Spdlines_sltol_%.1f_jp_%d_%s.vtk',cfg.out_dir,outdir,tag,glacfg.sl_tol,glacfg.jp,outversion),'polydata','lines',SL_tmp(1:glacfg.jp:end),S_tmp(1:glacfg.jp:end));
vtkwrite_spdlines(sprintf('%s/%s/%s_Rlines_sltol_%.1f_jp_%d_%s.vtk',cfg.out_dir,outdir,tag,glacfg.sl_tol,glacfg.jp,outversion),'polydata','lines',SL_tmp(1:glacfg.jp:end),R_tmp(1:glacfg.jp:end));
vtkwrite(sprintf('%s/%s/%s_disp_lentol_%.2f_%s.vtk',cfg.out_dir,outdir,tag,glacfg.sl_tol,outversion), 'structured_grid', PATH.startp(indb,1), PATH.startp(indb,2), PATH.startp(indb,3), ... 
    'vectors', 'vector_field', PATH.disp(indb,1), PATH.disp(indb,2), PATH.disp(indb,3));
if isfield(cfg,'anato')
vtkwrite(sprintf('%s/%s/%s_anato_%s.vtk',cfg.out_dir,outdir,tag,outversion), 'structured_points', 'mask', anato);
end

%%
T = toc;
fprintf('\n Post Processing --- Elapsed Time: %s\n',datestr(seconds(T),'HH:MM:SS'))

%%

Lag.s = s;
Lag.Pe = Pe;
Lag.r = r;
Lag.s_full = s_full;
Lag.Pe_rull = Pe_full;
Lag.r_rull = r_full;

Lag.PATH = PATH;
Lag.SL = SL;
Lag.sstream = sstream2;
Lag.Rstream = Rstream2;
Lag.pestream = pestream2;

end