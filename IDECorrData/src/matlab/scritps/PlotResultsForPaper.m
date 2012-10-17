% plot results for paper
close all

load Results_LinearModel_PBC_Aniso

%%
height = 6;
width = 18.2;
fig = figure('units','centimeters','position',[0 0 width height],'filename',filename,...
    'papersize',[height, width],'paperorientation','landscape','renderer','painters');

ax1 = axes('parent',fig,'units','centimeters','position',[1.2 1.5 3 3]);

imagesc(r,r,w,[cmin,cmax])
xlabel('Space','fontsize',FS_Label)
ylabel('Space','fontsize',FS_Label)
xlim([-10,10])
ylim([-10,10])
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10],'fontsize',FS_Tick)
axis square
axis xy
colorbar('units','centimeters','location','northoutside','position',[1.2 4.8 3 0.5],'fontsize',FS_Label)
text(-1, -18,'(a)','fontsize',FS_Label,'fontname','times roman')


spacing = 4.4;
ax2 = axes('parent',fig,'units','centimeters','position',[1.2+spacing 1.5 3 3]);

r_est = linspace(-20,20,size(w_est,1));
imagesc(r_est,r_est,(w_est),[cmin,cmax])
xlabel('Space','fontsize',FS_Label)
% ylabel('Space','fontsize',FS_Label)
xlim([-10,10])
ylim([-10,10])
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10],'fontsize',FS_Tick)
axis square
axis xy
colorbar('units','centimeters','location','northoutside','position',[1.2+spacing 4.8 3 0.5],'fontsize',FS_Label)
text(-1, -18,'(b)','fontsize',FS_Label,'fontname','times roman')



ax4 = axes('parent',fig,'units','centimeters','position',[1.2+3*spacing 1.5 3 3]);

plot(r,w(21,:),'k'), hold on
plot(r_est,w_est(14,:),'+r')


xlim([-10,10])
ylim([cmin,cmax])
xlabel('Space','fontsize',FS_Label)
ylabel('Kernel Amplitude','fontsize',FS_Label)
box off
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10.0],'fontsize',FS_Tick)
text(-1, -22,'(d)','fontsize',FS_Label,'fontname','times roman')


load Results_NonlinearModel_PBC_Aniso
plot(r_est,w_est(14,:),'xb')
hold off

leg = legend('$w(\mathbf{r})$','$\hat{w}_{L} (\mathbf{r})$', '$\hat{w}_{N} (\mathbf{r})$');
set(leg,'interpreter','latex','box','off','position',[1.2+3*spacing 5 5 0.15],'orientation','vertical','fontsize',FS_Label)

ax3 = axes('parent',fig,'units','centimeters','position',[1.2+2*spacing 1.5 3 3]);

r_est = linspace(-20,20,size(w_est,1));
imagesc(r_est,r_est,(w_est),[cmin,cmax])
xlabel('Space','fontsize',FS_Label)
% ylabel('Space','fontsize',FS_Label)
xlim([-10,10])
ylim([-10,10])
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10],'fontsize',FS_Tick)
axis square
axis xy
colorbar('units','centimeters','location','northoutside','position',[1.2+2*spacing 4.8 3 0.5],'fontsize',FS_Label)
text(-1, -18,'(c)','fontsize',FS_Label,'fontname','times roman')


%%
clear
load Results_LinearModel_PBC_Iso

height = 6;
width = 18.2;
fig = figure('units','centimeters','position',[0 0 width height],'filename',filename,...
    'papersize',[height, width],'paperorientation','landscape','renderer','painters');

ax1 = axes('parent',fig,'units','centimeters','position',[1.2 1.5 3 3]);

imagesc(r,r,w,[cmin,cmax])
xlabel('Space','fontsize',FS_Label)
ylabel('Space','fontsize',FS_Label)
xlim([-10,10])
ylim([-10,10])
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10],'fontsize',FS_Tick)
axis square
axis xy
colorbar('units','centimeters','location','northoutside','position',[1.2 4.8 3 0.5],'fontsize',FS_Label)
text(-1, -18,'(a)','fontsize',FS_Label,'fontname','times roman')


spacing = 4.4;
ax2 = axes('parent',fig,'units','centimeters','position',[1.2+spacing 1.5 3 3]);

r_est = linspace(-20,20,size(w_est,1));
imagesc(r_est,r_est,(w_est),[cmin,cmax])
xlabel('Space','fontsize',FS_Label)
% ylabel('Space','fontsize',FS_Label)
xlim([-10,10])
ylim([-10,10])
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10],'fontsize',FS_Tick)
axis square
axis xy
colorbar('units','centimeters','location','northoutside','position',[1.2+spacing 4.8 3 0.5],'fontsize',FS_Label)
text(-1, -18,'(b)','fontsize',FS_Label,'fontname','times roman')



ax4 = axes('parent',fig,'units','centimeters','position',[1.2+3*spacing 1.5 3 3]);

plot(r,w(21,:),'k'), hold on
plot(r_est,w_est(14,:),'+r')


xlim([-10,10])
ylim([-10,27])
xlabel('Space','fontsize',FS_Label)
ylabel('Kernel Amplitude','fontsize',FS_Label)
box off
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 25.0],'fontsize',FS_Tick)
text(-1, -22,'(d)','fontsize',FS_Label,'fontname','times roman')


load Results_NonlinearModel_PBC_Iso
plot(r_est,w_est(14,:),'xb')
hold off

leg = legend('$w(\mathbf{r})$','$\hat{w}_{L} (\mathbf{r})$', '$\hat{w}_{N} (\mathbf{r})$');
set(leg,'interpreter','latex','box','off','position',[1.2+3*spacing 5 5 0.15],'orientation','vertical','fontsize',FS_Label)

ax3 = axes('parent',fig,'units','centimeters','position',[1.2+2*spacing 1.5 3 3]);

r_est = linspace(-20,20,size(w_est,1));
imagesc(r_est,r_est,(w_est),[cmin,cmax])
xlabel('Space','fontsize',FS_Label)
% ylabel('Space','fontsize',FS_Label)
xlim([-10,10])
ylim([-10,10])
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10],'fontsize',FS_Tick)
axis square
axis xy
colorbar('units','centimeters','location','northoutside','position',[1.2+2*spacing 4.8 3 0.5],'fontsize',FS_Label)
text(-1, -18,'(c)','fontsize',FS_Label,'fontname','times roman')
