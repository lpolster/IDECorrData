%# create stacked images (I am simply repeating the same image 5 times)
% img = load('clown');
% I = repmat(img.X,[1 1 5]);
% cmap = img.map;
% 
% imagesc(cmap)
clc
close all
clear

load Results_NonlinearModel_PBC_Aniso

%%

FS = 8;
%# coordinates
[X,Y] = meshgrid(1:41, 1:41);
Z = ones(41:41);

height = 6;
width = 8.8;
fig = figure('units','centimeters','position',[0 0 width height],'filename',filename,...
    'papersize',[height, width],'paperorientation','landscape','renderer','painters');
subplot(121)
%# plot each slice as a texture-mapped surface (stacked along the Z-dimension)
for k=100:2:110
    surface('XData',X-0.5, 'YData',Y-0.5, 'ZData',Z.*k, ...
        'CData',squeeze(v(k,:,:)),'linestyle','none')
    shading interp
%     
%     , 'CDataMapping','direct', ...
%         'EdgeColor','none', 'FaceColor','texturemap')
end
% colormap(hot)
view(3), box off, axis tight
% set(gca, 'YDir','reverse', 'ZLim',[0 size(I,3)+1])
zlabel('Time','fontsize',FS)
ylabel('Space','fontsize',FS)
xlabel('Space','fontsize',FS)
% set(gca,'fontsize',FS,'xtick',[1,21,41],'xticklabel',{'-10','0','10'},'ytick',[1,21,41],'yticklabel',{'-10','0','10'})
set(gca,'fontsize',FS,'xtick',[],'ytick',[],'ztick',[])
view([-45 15])

subplot(324)
surf(w,'linestyle','none','facealpha',0.5)
axis tight
shading interp
set(gca,'fontsize',FS,'xtick',[],'ytick',[],'ztick',[])
view([-45 15])
zlabel('Kernel Amplitude','fontsize',FS)
ylabel('Space','fontsize',FS)
xlabel('Space','fontsize',FS)

%%
clc
close all

FS = 8;
%# coordinates
[X,Y] = meshgrid(1:41, 1:41);
Z = ones(41:41);

height = 10;
width = 8.8;

fig = figure('units','centimeters','position',[0 0 width height],...
    'papersize',[height, width],'paperorientation','landscape','renderer','painters');

for k=200:300
    surface('XData',X-0.5, 'YData',Y-0.5, 'ZData',Z.*k, ...
        'CData',squeeze(v(k,:,:)),'linestyle','none')
    shading interp
end
axis tight
box on
% set(gca, 'YDir','reverse', 'ZLim',[0 size(I,3)+1])
zlabel('Time (ms)','fontsize',FS)
ylabel('Space','fontsize',FS)
xlabel('Space','fontsize',FS)
% set(gca,'fontsize',FS,'xtick',[1,21,41],'xticklabel',{'-10','0','10'},'ytick',[1,21,41],'yticklabel',{'-10','0','10'})
set(gca,'fontsize',FS,'xtick',[1 21 41],'xticklabel',{'-10','0','10'},'ytick',[1 21 41],'yticklabel',{'-10','0','10'},'ztick',[1 50 100],'zticklabel',{'0','50','100'})
view([-45 45])

%%

height = 6;
width = 8.8;

fig = figure('units','centimeters','position',[0 0 width height],'filename',filename,...
    'papersize',[height, width],'paperorientation','landscape','renderer','painters');
