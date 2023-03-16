% This script shows how to run urOMT algorithm and Eulerian & Lagrangian post-processing on a rat brain dataset
% Parameters of urOMT are defined in getPamams.m
% Parameters of Eulerian post-processing are defined in paramInitEULApar.m
% Parameters of Lagrangian post-processing are defined in paramInitGLADpar.m
% One can run this script directly with default parameters.
% This script comes in 5 parts (see below)
% Created by Xinan Chen in March 2023

%%
addpath('utilities', 'Sensitivities','Inverse');

tag = 'C1217';


cfg = getParams(tag);

%% Part 1: load in data & mask    
if cfg.mask_number
    tmp = load_rhon(cfg.ROI_msk_path);
    tmp = tmp(cfg.x_range,cfg.y_range,cfg.z_range);
    cfg.msk = zeros(size(tmp));
    cfg.msk(tmp>0) = 1; 
    cfg.msk = bwmorph3(cfg.msk,'fill');
else
    cfg.msk = ones(cfg.true_size);
end
if cfg.do_resize
   cfg.msk = resizeMatrix(double(cfg.msk),round(cfg.size_factor.*size(cfg.msk)),'linear');
   cfg.msk(cfg.msk~=1) = 0;
end

% load data and chi if any
fprintf('--- Loading in data...\n')
if isfield(cfg,'chi_dir')
    imgs = getData(tag,cfg.first_time:cfg.time_jump:cfg.last_time+cfg.time_jump,cfg.first_time:cfg.time_jump:cfg.last_time);
else
    imgs = getData(tag,cfg.first_time:cfg.time_jump:cfg.last_time+cfg.time_jump);
end
for j = 1:(cfg.last_time-cfg.first_time)/cfg.time_jump+2
    tmp = imgs.im_data(:,:,:,j);
    if cfg.smooth>0
       tmp = affine_diffusion_3d(tmp,cfg.smooth,0.1,1,1);
    end
    if cfg.do_resize
       tmp = resizeMatrix(double(tmp),round(cfg.size_factor.*size(tmp)),'linear');
    end
    tmp(cfg.msk==0) = 0;
    cfg.vol(j).data = tmp;
end
if isfield(cfg,'chi_dir')
for j = 1:(cfg.last_time-cfg.first_time)/cfg.time_jump+1
    tmp = imgs.im_chi(:,:,:,:,j);
    if cfg.do_resize
       tmp = resizeMatrix(double(tmp),round(cfg.size_factor.*size(tmp)),'linear');
    end
    tmp(cfg.msk==0) = 0;
    cfg.chi(j).data = tmp;
end
end
fprintf('--- Data loaded.\n\n')


%% Part 2: Run urOMT

[cfg, flag] = runUROMT(cfg);

%% Part 3: Run Eulerian post-processing

[cfg, Eul] = runEULA(cfg);


%% Part 4: Run Lagrangian post-processing

[cfg, Lag] = runGLAD(cfg);


%% Part 5: Visualizations


%% Eul

% Eul averaged speed
x = 1:cfg.true_size(1);
y = 1:cfg.true_size(2);
z = 1:cfg.true_size(3);

for i = 1:2:length(cfg.first_time:cfg.time_jump:cfg.last_time)
figure,
tmp = reshape(mean(Eul.Es(:,(i-1)*cfg.nt+1:(i+1)*cfg.nt),2),cfg.true_size);
hs=slice(y,x,z,tmp,y,x,z); 
set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
alpha('color'),alphamap(linspace(0,1,100))
grid off, box off, axis image
xticks(0:10:cfg.true_size(2)); yticks(0:10:cfg.true_size(1)); zticks(0:10:cfg.true_size(3));
ax = gca; ax.FontSize = 10; 
%xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
colormap(jet)
clim([0,0.3])
view([-180   0])
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
xlim([0 cfg.true_size(2)]); ylim([0 cfg.true_size(1)]); zlim([0 cfg.true_size(3)])

saveas(gcf, sprintf('%s/%s/%s_EulAveSpeed_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Eul,cfg.tag,cfg.first_time+(i-1)*cfg.time_jump,cfg.first_time+(i+1)*cfg.time_jump)); 

end

% Eul averaged r
for i = 1:2:length(cfg.first_time:cfg.time_jump:cfg.last_time)
figure,
tmp = reshape(mean(Eul.Er(:,(i-1)*cfg.nt+1:(i+1)*cfg.nt),2),cfg.true_size);
tmp(tmp==0) = NaN;
hs=slice(y,x,z,tmp,y,x,z); 
set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
alpha('color'),alphamap([linspace(0,0.2,100)])
grid off, box off, axis image
xticks(0:10:cfg.true_size(2)); yticks(0:10:cfg.true_size(1)); zticks(0:10:cfg.true_size(3));
ax = gca; ax.FontSize = 10; 
%xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
clim([-0.04,0.04])
colormap(bluewhitered_drdt)
view([-180   0])
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
xlim([0 cfg.true_size(2)]); ylim([0 cfg.true_size(1)]); zlim([0 cfg.true_size(3)])

saveas(gcf, sprintf('%s/%s/%s_EulAveR_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Eul,cfg.tag,cfg.first_time+(i-1)*cfg.time_jump,cfg.first_time+(i+1)*cfg.time_jump)); 

end
%% Lag

%
% Lag pathlines

figure,
SL2 = Lag.SL(Lag.PATH.ind_msk);
nSL = length(SL2);
%colors = jet(round(max(PATH.displen)));

for ind = 1:10:nSL
    SL_tmp = SL2{ind};
    colors = jet(size(SL_tmp,1));
    hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
    %set(hlines,'EdgeColor',colors(round(PATH.displen(ind)),:));
    set(hlines,'FaceVertexCData',[colors;colors(end,:)],'EdgeColor','flat','FaceColor','none');
end
%
hold on;
axis image
colormap('jet');
view([-180   0])
xticks(0:10:cfg.true_size(2)); yticks(0:10:cfg.true_size(1)); zticks(0:10:cfg.true_size(3));
ax = gca; ax.FontSize = 10; 
%xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.3201 0.2515 0.4894 0.3880],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
grid off,
xlim([0 cfg.true_size(2)]); ylim([0 cfg.true_size(1)]); zlim([0 cfg.true_size(3)])

cb = colorbar;
set(cb,'YTick',[])
text(1,-3,53,'end','Rotation',90,'FontSize',25);
text(1,-3,-1,'start','Rotation',90,'FontSize',25);

saveas(gcf, sprintf('%s/%s/%s_LagPathlines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Lag_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 

%
% Lag Speed-lines
figure,
spdmax = 0.3;
Num = 100; % the larger, the more intervals in colors

SL2 = Lag.SL(Lag.PATH.ind_msk);
SL_spd = Lag.sstream(Lag.PATH.ind_msk);
nSL = length(SL2);
colors = jet(Num);

for ind = 1:10:nSL
    SL_tmp = SL2{ind};
    SL_spd_tmp = SL_spd{ind};
    SL_spd_tmp(SL_spd_tmp>spdmax) = spdmax;
    SL_spd_rk = round(SL_spd_tmp/spdmax*Num);
    SL_spd_rk(SL_spd_rk<1) = 1;
    hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
    set(hlines,'FaceVertexCData',[colors(SL_spd_rk,:);colors(SL_spd_rk(end),:)],'EdgeColor','flat','FaceColor','none');
end
%
hold on;
view([-180   0])
xticks(0:10:cfg.true_size(2)); yticks(0:10:cfg.true_size(1)); zticks(0:10:cfg.true_size(3));
ax = gca; ax.FontSize = 10; 
%xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.3201 0.2515 0.4894 0.3880],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
grid off; axis image;
xlim([0 cfg.true_size(2)]); ylim([0 cfg.true_size(1)]); zlim([0 cfg.true_size(3)])
colormap('jet'); grid off;

cb = colorbar;
cb.Ticks = linspace(0, 1, 4);
cb.TickLabels = num2cell(0:0.1:0.3);
cb.FontSize = 20;
text(1,-3,15,'speed (a.u.)','Rotation',90,'FontSize',25);

xlim([0 cfg.true_size(2)]); ylim([0 cfg.true_size(1)]); zlim([0 cfg.true_size(3)])
saveas(gcf, sprintf('%s/%s/%s_LagSpdlines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Lag_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 


%
% Lag Peclet-lines
figure,
Pemax = 250;
Num = 100; % the larger, the more intervals in colors

SL2 = Lag.SL(Lag.PATH.ind_msk);
SL_Pe = Lag.pestream(Lag.PATH.ind_msk);
nSL = length(SL2);
colors = jet(Num);

for ind = 1:10:nSL
    SL_tmp = SL2{ind};
    SL_Pe_tmp = SL_Pe{ind};
    SL_Pe_tmp(SL_Pe_tmp>Pemax) = Pemax;
    SL_Pe_rk = round(SL_Pe_tmp/Pemax*Num);
    SL_Pe_rk(SL_Pe_rk<1) = 1;
    hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
    set(hlines,'FaceVertexCData',[colors(SL_Pe_rk,:);colors(SL_Pe_rk(end),:)],'EdgeColor','flat','FaceColor','none');
end
%
hold on;
view([-180   0])
xticks(0:10:cfg.true_size(2)); yticks(0:10:cfg.true_size(1)); zticks(0:10:cfg.true_size(3));
ax = gca; ax.FontSize = 10; 
%xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.3201 0.2515 0.4894 0.3880],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
grid off; axis image;
xlim([0 cfg.true_size(2)]); ylim([0 cfg.true_size(1)]); zlim([0 cfg.true_size(3)])
colormap('jet'); grid off;

cb = colorbar;
cb.Ticks = linspace(0, 1, 6);
cb.TickLabels = num2cell(0:50:250);
cb.FontSize = 20;
text(1,-3,25,'\it{Pe}','Rotation',90,'FontSize',25);

xlim([0 cfg.true_size(2)]); ylim([0 cfg.true_size(1)]); zlim([0 cfg.true_size(3)])
saveas(gcf, sprintf('%s/%s/%s_LagPelines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Lag_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 

%
% Lag Velocity flux vectors
figure,
strid = 10;
magnify = 1;
[x, y, z] = meshgrid(1:cfg.true_size(2), 1:cfg.true_size(1), 1:cfg.true_size(3));
mskfv = isosurface(x,y,z,cfg.msk,0.5);
mskp = patch(mskfv);
mskp.FaceColor = [.17,.17,.17];
mskp.FaceAlpha= 0.031;
mskp.EdgeColor = [.17,.17,.17];
mskp.EdgeAlpha= 0;
mskp.DisplayName = 'mask';

grid off, axis image
hold on,
q = quiver3(Lag.PATH.startp(Lag.PATH.ind_msk(1:strid:end),2),Lag.PATH.startp(Lag.PATH.ind_msk(1:strid:end),1),Lag.PATH.startp(Lag.PATH.ind_msk(1:strid:end),3),Lag.PATH.disp(Lag.PATH.ind_msk(1:strid:end),2)*magnify,Lag.PATH.disp(Lag.PATH.ind_msk(1:strid:end),1)*magnify,Lag.PATH.disp(Lag.PATH.ind_msk(1:strid:end),3)*magnify,...
    'color','r','LineWidth',1.5,'MaxHeadSize',0.5,'AutoScale','off','DisplayName','flux vectors');
%title(sprintf('Velocity Flux Vectors'),'FontSize',20, 'Interpreter', 'none'), 
view([-180   0])
ax = gca; ax.FontSize = 10; 
%xlabel('x-axis','FontSize',18),ylabel('y-axis','FontSize',18),zlabel('z-axis','FontSize',18)
set(get(gca,'YLabel'),'Rotation',10);%xtickangle(20);
set(get(gca,'XLabel'),'Rotation',-30,'VerticalAlignment','middle');
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.3201 0.2515 0.4894 0.3880],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
grid off; axis image;
xticks(0:10:cfg.true_size(2)); yticks(0:10:cfg.true_size(1)); zticks(0:10:cfg.true_size(3));
%// Compute the magnitude of the vectors
mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
            reshape(q.WData, numel(q.UData), [])).^2, 2));

%// Get the current colormap
currentColormap = colormap(jet);
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   
set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');

view([-180   0])
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.3201 0.2515 0.4894 0.3880],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
%colorbar('FontSize',18), 
grid off,
xlim([0 cfg.true_size(2)]); ylim([0 cfg.true_size(1)]); zlim([0 cfg.true_size(3)])

cb = colorbar;
cb.Ticks = linspace(0, 800, 6);
cb.TickLabels = num2cell(0:160:800);
cb.FontSize = 20;
text(1,-4,-10,'distance of movement (a.u.)','Rotation',90,'FontSize',25);

saveas(gcf, sprintf('%s/%s/%s_LagFluxVector_E%02d_%02d.png',cfg.out_dir,cfg.outdir_Lag_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 


% input images
x = 1:cfg.true_size(1);
y = 1:cfg.true_size(2);
z = 1:cfg.true_size(3);
for i = cfg.first_time:cfg.time_jump:cfg.last_time+cfg.time_jump
    figure,
    tmp = cfg.vol((i-cfg.first_time)/cfg.time_jump+1).data;
    hs=slice(y,x,z,tmp,y,x,z); 
    set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
    alpha('color'),alphamap(gca,linspace(0,1,100))
    xticks(0:10:cfg.true_size(2)); yticks(0:10:cfg.true_size(1)); zticks(0:10:cfg.true_size(3));
    ax = gca; ax.FontSize = 10; 
    clim([0,300])
    colormap(gca,jet)
    view([-180   0])
    set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5152 0.2790 0.2943 0.3605],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
    grid off, axis image
    xlim([0 cfg.true_size(2)]); ylim([0 cfg.true_size(1)]); zlim([0 cfg.true_size(3)])

    saveas(gcf,sprintf('%s/images_%d_tind_%d.png',cfg.out_dir,i,(i-cfg.first_time)/cfg.time_jump+1))
end



