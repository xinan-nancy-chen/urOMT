function glacfg = paramInitGLADpar(cfg)
% Created by Xinan Chen on 04/26/2021.
% This function returns parameters corresponding to the input 'data_name' and analysis type
%   Output: 'gladcfg' := structure containing parameters for GLAD analysis

%%

glacfg.do_masked = 1;%1 if ensure everything is inside mask before saving to .nii
glacfg.spMSK_str = 'ROImsk'; % original ROI mask for rOMT

%% set parameters, start points and directory
switch cfg.dataset_name
    case {'RatBrainsCAA3M'}
        glacfg.do_sp = 1; %1 if use sp from max(data)>sp_thresh, 0 o/w
        glacfg.sp_thresh = 10;%6;%8;%1;%12;

        glacfg.sl_tol = 1;%2; %threshold for minimum Euclidean length between initial and final streamline points                
        glacfg.spType = 'ordered'; %'uniform' (default: 'ordered'); %'ordered': every fs-th possible sp is used (default)sp distribution; 'uniform': sp chosen according to uniform distribution
        switch glacfg.spType
            case 'ordered'
                glacfg.spfs = 1;%5; % start point jump
                glacfg.spSTR = sprintf('spDISTordered_fs%d',glacfg.spfs);
                glacfg.spSTRtitle = sprintf('spDISTordered-fs%d',glacfg.spfs);
            case 'uniform'
                glacfg.spPerc = 40;
                glacfg.spSTR = sprintf('spDISTunif_spPerc%d',glacfg.spPerc);
                glacfg.spSTRtitle = sprintf('spDISTunif-spPerc%d',glacfg.spPerc);
        end
        glacfg.nEulStep = 1; %number of Eulerian steps to be taken
        
        glacfg.smoothv = 0; % smooth velocity field
        if glacfg.smoothv
            glacfg.Svt = 5;%10; % smoothing w.r.t time
            glacfg.Svs = 1;%5; % smoothing w.r.t space
            glacfg.smoothvSTR = sprintf('1_Tolt%d_Tols%d',glacfg.Svt,glacfg.Svs);
        else
            glacfg.smoothvSTR = '0';
        end
        
        glacfg.smoothp = 0;%1; % smooth pathline
        if glacfg.smoothp
            glacfg.smpTol = 100; % the higher, the smoother the pathlines
            glacfg.smoothpSTR = sprintf('1_Tol%d',glacfg.smpTol);
        else
            glacfg.smoothpSTR = '0';
        end
        
        glacfg.do_QB = 0;%1;%0
        if glacfg.do_QB
            glacfg.mdt = 10;
            glacfg.sfs = 5;%1;%5;
            % QuickBundle parameters 
            glacfg.qbthresh = 4;%5;%QuickBundles threshold distance
            glacfg.clus_tol = 12; %threshold for min # of sl's needed in a cluster
            glacfg.clus_cutoff = 600; %pick up to # of largest clusters
            glacfg.nbp = 124; %# of points that QuickBundles will resample streamline to
            glacfg.metric_str = 'AveragePointwiseEuclideanMetric';
            switch glacfg.metric_str
                case 'AveragePointwiseEuclideanMetric'
                    glacfg.ms = 'AvgPwEuc';
                    glacfg.feature_str = sprintf('ResampleFeature(nb_points=%d)',glacfg.nbp);
                case 'CosineMetric'
                    glacfg.ms = 'Cos';
                    glacfg.feature_str = 'VectorOfEndpointsFeature()';
            end
            glacfg.QB_str = sprintf('QB_qbthresh%d_clus_tol%d_clus_cutoff%d_%s',glacfg.qbthresh,glacfg.clus_tol,glacfg.clus_cutoff,glacfg.feature_str);
        else
            glacfg.QB_str = 'QB0';
            glacfg.mdt = 1;
            glacfg.sfs = 1;%5;
        end
        
        glacfg.vfs = 1;

        % vis parameters
        glacfg.jp = 5;%1;%5;

        glacfg.cutoff_str = 'min';
        switch glacfg.cutoff_str
            case 'min'
                glacfg.thresholds.conc = 0.01;%0.0001;
                glacfg.thresholds.front = 0;
                glacfg.thresholds.flow = 0;
                glacfg.thresholds.speed = 0.01;%0.0001;
            case 'max'
                glacfg.thresholds.conc = .0001;
                glacfg.thresholds.front = 0;
                glacfg.thresholds.flow = 0;
                glacfg.thresholds.speed = 0;
            case 'mean'
                glacfg.thresholds.conc = .00001;
                glacfg.thresholds.front = 0;
                glacfg.thresholds.flow = 0.0001;
                glacfg.thresholds.speed = 0.0001;
        end

    case {'Gauss'}
        glacfg.do_sp = 1; %1 if use sp from max(data)>sp_thresh, 0 o/w
        glacfg.sp_thresh = 5;%6;%8;%1;%12;

        glacfg.sl_tol = 1;%2; %threshold for minimum Euclidean length between initial and final streamline points                
        glacfg.spType = 'ordered'; %'uniform' (default: 'ordered'); %'ordered': every fs-th possible sp is used (default)sp distribution; 'uniform': sp chosen according to uniform distribution
        switch glacfg.spType
            case 'ordered'
                glacfg.spfs = 1;%5; % start point jump
                glacfg.spSTR = sprintf('spDISTordered_fs%d',glacfg.spfs);
                glacfg.spSTRtitle = sprintf('spDISTordered-fs%d',glacfg.spfs);
            case 'uniform'
                glacfg.spPerc = 40;
                glacfg.spSTR = sprintf('spDISTunif_spPerc%d',glacfg.spPerc);
                glacfg.spSTRtitle = sprintf('spDISTunif-spPerc%d',glacfg.spPerc);
        end
        glacfg.nEulStep = 1; %number of Eulerian steps to be taken
        
        glacfg.smoothv = 0; % smooth velocity field
        if glacfg.smoothv
            glacfg.Svt = 5;%10; % smoothing w.r.t time
            glacfg.Svs = 1;%5; % smoothing w.r.t space
            glacfg.smoothvSTR = sprintf('1_Tolt%d_Tols%d',glacfg.Svt,glacfg.Svs);
        else
            glacfg.smoothvSTR = '0';
        end
        
        glacfg.smoothp = 0;%1; % smooth pathline
        if glacfg.smoothp
            glacfg.smpTol = 100; % the higher, the smoother the pathlines
            glacfg.smoothpSTR = sprintf('1_Tol%d',glacfg.smpTol);
        else
            glacfg.smoothpSTR = '0';
        end
        
        glacfg.do_QB = 0;%1;%0
        if glacfg.do_QB
            glacfg.mdt = 10;
            glacfg.sfs = 5;%1;%5;
            % QuickBundle parameters 
            glacfg.qbthresh = 4;%5;%QuickBundles threshold distance
            glacfg.clus_tol = 12; %threshold for min # of sl's needed in a cluster
            glacfg.clus_cutoff = 600; %pick up to # of largest clusters
            glacfg.nbp = 124; %# of points that QuickBundles will resample streamline to
            glacfg.metric_str = 'AveragePointwiseEuclideanMetric';
            switch glacfg.metric_str
                case 'AveragePointwiseEuclideanMetric'
                    glacfg.ms = 'AvgPwEuc';
                    glacfg.feature_str = sprintf('ResampleFeature(nb_points=%d)',glacfg.nbp);
                case 'CosineMetric'
                    glacfg.ms = 'Cos';
                    glacfg.feature_str = 'VectorOfEndpointsFeature()';
            end
            glacfg.QB_str = sprintf('QB_qbthresh%d_clus_tol%d_clus_cutoff%d_%s',glacfg.qbthresh,glacfg.clus_tol,glacfg.clus_cutoff,glacfg.feature_str);
        else
            glacfg.QB_str = 'QB0';
            glacfg.mdt = 1;
            glacfg.sfs = 1;%5;
        end
        
        glacfg.vfs = 1;

        % vis parameters
        glacfg.jp = 5;%1;%5; 

        glacfg.cutoff_str = 'min';
        switch glacfg.cutoff_str
            case 'min'
                glacfg.thresholds.conc = 0.0001;
                glacfg.thresholds.front = 0;
                glacfg.thresholds.flow = 0;
                glacfg.thresholds.speed = 0.0001;
            case 'max'
                glacfg.thresholds.conc = .0001;
                glacfg.thresholds.front = 0;
                glacfg.thresholds.flow = 0;
                glacfg.thresholds.speed = 0;
            case 'mean'
                glacfg.thresholds.conc = .00001;
                glacfg.thresholds.front = 0;
                glacfg.thresholds.flow = 0.0001;
                glacfg.thresholds.speed = 0.0001;
        end

    otherwise
        fprintf('In paramInitGLADpar.m: non-applicable dataset_name!');      

end


glacfg.flw_type = 'vel';%'flw';
glacfg.pln = 2; %minimum number of unique points needed to be considered a pathline
glacfg.RD = 'R'; %'D' if use data; 'R' if use interpolated img
glacfg.XT = 'T'; %'X' if interp fixed spatial dist, 'T' if use time step;
glacfg.minIm0 = 0; %1 if set img = img - min(img), 0 otherwise



end