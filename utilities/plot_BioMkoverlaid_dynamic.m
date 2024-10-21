%%
addpath(genpath('utilites/'))
%% Suppose one has run Part 1 ~ part 4 to get cfg and Eul, then use the following for dynamic visualization

%%
% load in image data
clear imgs;
imgs.im = cfg.vol;
imgs.speed = Eul.Es;%reshape(Eul.Es,[cfg.true_size,length(ti:tj:tf)*nt]);
imgs.flux = Eul.Es.*Eul.Eintp(:,1:end-1);
imgs.Peclet = Eul.Epe; 
imgs.rflux = Eul.Er.*Eul.Eintp(:,1:end-1);
imgs.r = Eul.Er;

% set visualization parameters
clear vis;
vis.RhoMax = 80; % upper limit for input concentration images
vis.PeMax = 400; % upper limit for Peclet
vis.SMax = 0.6; % upper limit for speed
vis.Fluxmax = 6; % upper limit for flux
vis.RFluxmax = 4; % upper limit for rflux
vis.RMaxmin = 0.1; % upper limit for r

vis.hgt = round(cfg.true_size(3)/2)+1; % select the z-dimension slice for visualization

vis.titlestr = 'Eul';
vis.titleend = '';
vis.colorbarfontsz = 9;
vis.titleftsz = 12;
vis.ti = cfg.first_time;
vis.tj = cfg.time_jump;
vis.tf = cfg.last_time;
vis.nt = cfg.nt;
vis.n = cfg.true_size;
vis.dynamic_jp = 1;
vis.do_permute = 0;
%vis.msk_tumor = cfg.msk>0;



figure, plot_BioMk_dynamic(imgs,vis)
set(gcf,'position',[30,91,1408,690],'Color',[1,1,1]);




function  plot_BioMk_dynamic(imgs,vis,do_thresh)
% Created by Xinan Chen on 02/14/2023
% To create subplots of 5 biomarker plots, speed, flux, Peclet, rflux

%%

if nargin<3
    do_thresh = 0;
end

Nrow = length(fieldnames(imgs));


colorbarfontsz = vis.colorbarfontsz;
RhoMax = vis.RhoMax;
RMaxmin = vis.RMaxmin;
PeMax = vis.PeMax;
SMax = vis.SMax;
RFluxmax = vis.RFluxmax;
Fluxmax = vis.Fluxmax;

titleftsz = vis.titleftsz;
ti = vis.ti;
tj = vis.tj;
tf = vis.tf;
nt = vis.nt;
n = vis.n;
hgt = vis.hgt;
jp = vis.dynamic_jp;

range = ti:tj:tf+tj;
range = range(1:jp:end);
order = 1:length(ti:tj:tf+tj);
order = order(1:jp:end);


EulsAll = zeros([n,length(1:jp:length(ti:tj:tf))]);
EulfluxAll = EulsAll;
EulPeAll = EulsAll;
EulrAll = EulsAll;
EulrfluxAll = EulsAll;

for i = 1:jp:length(ti:tj:tf)
    EulsAll(:,:,:,i) = reshape(sum(imgs.speed(:,(i-1)*nt+1:(i-1+jp)*nt),2),n)/(jp*nt);
    EulfluxAll(:,:,:,i) = reshape(sum(imgs.flux(:,(i-1)*nt+1:(i-1+jp)*nt),2),n)/(jp*nt);
    EulPeAll(:,:,:,i) = reshape(sum(imgs.Peclet(:,(i-1)*nt+1:(i-1+jp)*nt),2),n)/(jp*nt);
    EulrAll(:,:,:,i) = reshape(sum(imgs.r(:,(i-1)*nt+1:(i-1+jp)*nt),2),n)/(jp*nt);
    EulrfluxAll(:,:,:,i) = reshape(sum(imgs.rflux(:,(i-1)*nt+1:(i-1+jp)*nt),2),n)/(jp*nt);
end

im = imgs.im;
speed = EulsAll;
flux = EulfluxAll;
Peclet = EulPeAll;
r = EulrAll;
rflux = EulrfluxAll;

%% 

NN = 1;
% input
for i = range
    subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01,'Margin',0.04),
    tmp = im(NN).data;
    %tmp(msk_tumor==0) = 0;
    if vis.do_permute
        tmp = permute(tmp,[2,1,3]);
    end
    montageArray(tmp(:,:,hgt));
    clim([0 RhoMax]);
    colormap(gca,gray(256));%colormap(gca,rgb2gray(colormap));
    axis off;
    axis image;
    title(sprintf('input\nF%d',i),'Fontsize',titleftsz,'Interpreter','none');
    NN = NN + 1;
end
subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01),
clim([0 RhoMax]);
axis off;
colormap(gca,gray(256));%colormap(gca,rgb2gray(colormap));
colorbar('Fontsize',colorbarfontsz)
hold on
NN = NN + 1;

% speed
for i = 1:length(range)-1
    subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01,'Margin',0.04),
    tmp = mean(speed(:,:,:,order(i):order(i+1)-1),4);%reshape(sum(speed(:,(i-1)*nt+1:(i-1+EulerianJump)*nt),2),n)/(EulerianJump*nt);
    %tmp(msk_tumor==0) = 0;
    if do_thresh
        img = (im(order(i)).data+im(order(i)+1).data)/2;
        tmp(img<do_thresh) = 0;
    end
    if vis.do_permute
        tmp = permute(tmp,[2,1,3]);
    end
    montageArray(tmp(:,:,hgt));
    title(sprintf('Speed\nF%d -> %d',range(i),range(i+1)),'Fontsize',12,'Interpreter','none')
    clim([0,SMax])
    colormap(gca,'turbo');
    axis off;
    axis image;
    NN = NN + 1;
end
NN = NN + 1;
subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01),
clim([0,SMax]);
axis off;
colormap(gca,'turbo');
colorbar('Fontsize',colorbarfontsz)
hold on
NN = NN + 1;

% flux
for i = 1:length(range)-1
    subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01,'Margin',0.04),
    tmp = mean(flux(:,:,:,order(i):order(i+1)-1),4);%flux(:,:,:,i);%reshape(sum(flux(:,(i-1)*nt+1:(i-1+EulerianJump)*nt),2),n)/(EulerianJump*nt);
    %tmp(msk_tumor==0) = 0;
    if do_thresh
        img = (im(order(i)).data+im(order(i)+1).data)/2;
        tmp(img<do_thresh) = 0;
    end
    if vis.do_permute
        tmp = permute(tmp,[2,1,3]);
    end
    montageArray(tmp(:,:,hgt));
    title(sprintf('Flux\nF%d -> %d',range(i),range(i+1)),'Fontsize',12,'Interpreter','none')
    clim([0,Fluxmax])
    colormap(gca,'turbo');
    axis off;
    axis image;
    NN = NN + 1;
end
NN = NN + 1;
subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01),
clim([0,Fluxmax]);
axis off;
colormap(gca,'turbo');
colorbar('Fontsize',colorbarfontsz)
hold on
NN = NN + 1;

% Peclet
for i = 1:length(range)-1
    subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01,'Margin',0.04),
    tmp = mean(Peclet(:,:,:,order(i):order(i+1)-1),4);%Peclet(:,:,:,i);%reshape(sum(Peclet(:,(i-1)*nt+1:(i-1+EulerianJump)*nt),2),n)/(EulerianJump*nt);
    %tmp(msk_tumor==0) = 0;
    if do_thresh
        img = (im(order(i)).data+im(order(i)+1).data)/2;
        tmp(img<do_thresh) = 0;
    end
    if vis.do_permute
        tmp = permute(tmp,[2,1,3]);
    end
    montageArray(tmp(:,:,hgt));
    title(sprintf('Peclet\nF%d -> %d',range(i),range(i+1)),'Fontsize',12,'Interpreter','none')
    clim([0,PeMax])
    colormap(gca,'turbo');
    axis off;
    axis image;
    NN = NN + 1;
end
NN = NN + 1;
subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01),
clim([0,PeMax]);
axis off;
colormap(gca,'turbo');
colorbar('Fontsize',colorbarfontsz)
hold on
NN = NN + 1;


% rflux
for i = 1:length(range)-1
    subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01,'Margin',0.04),
    tmp = mean(rflux(:,:,:,order(i):order(i+1)-1),4);%rflux(:,:,:,i);%reshape(sum(rflux(:,(i-1)*nt+1:(i-1+EulerianJump)*nt),2),n)/(EulerianJump*nt);
    %tmp(msk_tumor==0) = 0;
    if do_thresh
        img = (im(order(i)).data+im(order(i)+1).data)/2;
        tmp(img<do_thresh) = 0;
    end
    if vis.do_permute
        tmp = permute(tmp,[2,1,3]);
    end
    montageArray(tmp(:,:,hgt));
    title(sprintf('influx/efflux\nF%d -> %d',range(i),range(i+1)),'Fontsize',12,'Interpreter','none')
    clim([-RFluxmax,RFluxmax])
    colormap(gca,bluewhitered_drdt);
    axis off;
    axis image;
    NN = NN + 1;
end
NN = NN + 1;
subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01),
clim([-RFluxmax,RFluxmax]);
axis off;
colormap(gca,bluewhitered_drdt);
colorbar('Fontsize',colorbarfontsz)
%hold on
%NN = NN + 1;

% % r
% for i = 1:length(range)-1
%     subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01,'Margin',0.04),
%     tmp = mean(r(:,:,:,order(i):order(i+1)-1),4);%r(:,:,:,i);%reshape(sum(r(:,(i-1)*nt+1:(i-1+EulerianJump)*nt),2),n)/(EulerianJump*nt);
%     %tmp(msk_tumor==0) = 0;
%     if do_thresh
%         img = (im(order(i)).data+im(order(i)+1).data)/2;
%         tmp(img<do_thresh) = 0;
%     end
%     if vis.do_permute
%         tmp = permute(tmp,[2,1,3]);
%     end
%     montageArray(tmp(:,:,hgt));
%     title(sprintf('r\nF%d -> %d',range(i),range(i+1)),'Fontsize',12,'Interpreter','none')
%     clim([-RMaxmin,RMaxmin])
%     colormap(gca,bluewhitered_drdt);
%     axis off;
%     axis image;
%     NN = NN + 1;
% end
% NN = NN + 1;
% subaxis(Nrow,length(range)+1,NN,'SpacingVert',0.05,'SpacingHoriz',0.01),
% clim([-RMaxmin,RMaxmin]);
% axis off;
% colormap(gca,bluewhitered_drdt);
% colorbar('Fontsize',colorbarfontsz)



end

