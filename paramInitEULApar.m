function eulcfg = paramInitEULApar(cfg)
% Created by Xinan Chen on 09/24/2022.
% This function returns parameters corresponding to the input 'data_name' and analysis type
%   Input:  type  := 's' if run for speed map; 'v' if run for flux vectors
%   Output: 'eulcfg' := structure containing parameters for EULA analysis

%%
%% set parameters, start points and directory
switch cfg.dataset_name
    case {'RatBrainsCAA3M','BreastTumor'}
        eulcfg.smoothv = 0; % smooth velocity field
        if eulcfg.smoothv
            eulcfg.Svt = 38; % smoothing w.r.t time
            eulcfg.Svt_type = 'gaussian';
            eulcfg.Svs = 0;%5; % smoothing w.r.t space
            eulcfg.smoothvSTR = sprintf('1_Tolt%s%d_Tols%d',eulcfg.Svt_type,eulcfg.Svt,eulcfg.Svs);
        else
            eulcfg.smoothvSTR = '0';
        end
        
        eulcfg.smoothr = 0; % smooth relative source r
        if eulcfg.smoothr
            eulcfg.Srt = 38; % smoothing w.r.t time
            eulcfg.Srt_type = 'gaussian';
            eulcfg.Srs = 0; % smoothing w.r.t space
            eulcfg.smoothrSTR = sprintf('1_Tolt%s%d_Tols%d',eulcfg.Srt_type,eulcfg.Srt,eulcfg.Srs);
        else
            eulcfg.smoothrSTR = '0';
        end
        %eulcfg.EulerianJump = 1;

    case {'Gauss'}
        eulcfg.smoothv = 0; % smooth velocity field
        if eulcfg.smoothv
            eulcfg.Svt = 38; % smoothing w.r.t time
            eulcfg.Svt_type = 'gaussian';
            eulcfg.Svs = 0;%5; % smoothing w.r.t space
            eulcfg.smoothvSTR = sprintf('1_Tolt%s%d_Tols%d',eulcfg.Svt_type,eulcfg.Svt,eulcfg.Svs);
        else
            eulcfg.smoothvSTR = '0';
        end
        
        eulcfg.smoothr = 0; % smooth relative source r
        if eulcfg.smoothr
            eulcfg.Srt = 38; % smoothing w.r.t time
            eulcfg.Srt_type = 'gaussian';
            eulcfg.Srs = 0; % smoothing w.r.t space
            eulcfg.smoothrSTR = sprintf('1_Tolt%s%d_Tols%d',eulcfg.Srt_type,eulcfg.Srt,eulcfg.Srs);
        else
            eulcfg.smoothrSTR = '0';
        end
        %eulcfg.EulerianJump = 1;

    otherwise
        fprintf('paramInitEULApar: non-applicable dataset_name!');     

end




end