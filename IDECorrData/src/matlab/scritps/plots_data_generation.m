% Dean Freestone. 
% test correlation analysis for estimation of kernel support

clc
clear
close all

UsePBC = true;
UseIso = true;
UseLinear = false;
UseLinearized = false;
UseNonlinear = true;


CheckDerivation = false;

tic                 % timer

% for plotting
% ~~~~~~~
FS_Label = 15;          % point / fontsize for the axis label 8
FS_Tick = 15;                % fontsize for the ticks 8
MS = 10;                     % marker size
LW = 1;
plotwidth_fig1 = 40;        % cm 4.4
plotheight_fig1 = 25;       %3.5
plotwidth_fig2 = 25;        
plotheight_fig2 = 15;
plotwidth_fig3 = 12;        
plotheight_fig3 = 12;
% parameters

%% spatial parameters
% ~~~~~~~~~~~
Delta = 0.5;                          % space step for the spatial discretisation
Delta_squared = Delta^2;
SpaceMax = 10;                    % maximum space in mm
SpaceMin = -SpaceMax;         % minimum space in mm
NPoints = (SpaceMax-SpaceMin)/Delta+1;
NPoints_total = NPoints^2;
r = linspace(SpaceMin,SpaceMax,NPoints);      % define space

EstimationSpaceMax = 10;
EstimationSpaceMin = -10;

%% temporal parameters
% ~~~~~~~~~~~~~~~~
Ts = 1e-3;              % sampling period (s)
T = 20000;              % maximum time (ms)          % 2 seconds = 100 seconds

%% spatial kernel parameters
% ~~~~~~~~~~~~~~

theta(1) = 100;%80.0;           % local kernel amplitude
theta(2) = -80;             % surround kernel amplitude
theta(3) = 5;              % lateral kernel amplitude
theta(4) = 15;              % anisotropic amplitude

sigma_psi(1) = 1.8;     % local kernel width
sigma_psi(2) = 2.4;     % surround kernel width
sigma_psi(3) = 6;       % lateral kernel width
sigma_psi(4) = 2;               % anisotropic width

psi_0 = Define2DGaussian(0,0, sigma_psi(1)^2, 0,NPoints,SpaceMin,SpaceMax);
psi_1 = Define2DGaussian(0,0, sigma_psi(2)^2, 0,NPoints,SpaceMin,SpaceMax);
psi_2 = Define2DGaussian(0,0, sigma_psi(3)^2, 0,NPoints,SpaceMin,SpaceMax);
psi_3 = Define2DGaussian(-3,0, sigma_psi(4)^2, 0,NPoints,SpaceMin,SpaceMax); % anisotropic

psi_0_scaled = theta(1)*psi_0;
psi_1_scaled = theta(2)*psi_1;
psi_2_scaled = theta(3)*psi_2;
psi_3_scaled = theta(4)*psi_3;

%w = theta(1)*psi_0 + theta(2)*psi_1 + theta(3)*psi_2+ theta(4)*psi_3;       % the kernel (k)
w = theta(1)*psi_0 + theta(2)*psi_1 + theta(3)*psi_2;  
%%
%plotting 3 scaled gaussians and resulting mexican hat
filename = 'C:\Documents and Settings\lpolster\IDECorrData\src\matlab\scritps\4KernelPlot.pdf';
figure('units','centimeters','position',[1 1 plotwidth_fig1 plotheight_fig1],'filename',filename,...
   'papersize',[plotheight_fig1, plotwidth_fig1],'paperorientation','landscape','renderer','painters') 
subplot(221)
surf(w)
title('Spatial Kernel','FontSize', 20)
xlabel('Space','FontSize', 15)
ylabel('Space','FontSize', 15)
zlabel('Connectivity Strength','FontSize', 15)
xlim([0,50])
ylim([0,50])
set(gca,'FontSize', 15)
colorbar('FontSize', 15)
subplot(222)
surf(psi_0_scaled)
title('Local Kernel','FontSize', 20)
xlabel('Space','FontSize', 15)
ylabel('Space','FontSize', 15)
zlabel('Connectivity Strength','FontSize', 15)
xlim([0,50])
ylim([0,50])
set(gca,'FontSize', 15)
colorbar('FontSize', 15)
subplot(223)
surf(psi_1_scaled)
title('Surround Kernel', 'FontSize', 20)
xlabel('Space','FontSize', 15)
ylabel('Space','FontSize', 15)
zlabel('Connectivity Strength','FontSize', 15)
xlim([0,50])
ylim([0,50])
set(gca,'FontSize', 15)
colorbar('FontSize', 15)
subplot(224)
surf(psi_2_scaled)
title('Lateral Kernel', 'FontSize', 20)
xlabel('Space','FontSize', 15)
ylabel('Space','FontSize', 15)
zlabel('Connectivity Strength','FontSize', 15)
xlim([0,50])
ylim([0,50])
set(gca,'xtick',[0 50],'ytick',[0 50],'FontSize', 15)
colorbar('FontSize', 15)
%%
%
x = linspace(SpaceMin,SpaceMax);  
psi_simple_0 = exp(-(((x).^2)./(sigma_psi(1)).^2));
psi_simple_1 = exp(-(((x).^2)./(sigma_psi(2)).^2));
psi_simple_2 = exp(-(((x).^2)./(sigma_psi(3)).^2));
psi_simple_3 = exp(-(((x+3).^2)./(sigma_psi(4)).^2));

psi_simple_0_scaled = theta(1)*psi_simple_0;
psi_simple_1_scaled = theta(2)*psi_simple_1;
psi_simple_2_scaled = theta(3)*psi_simple_2;
psi_simple_3_scaled = theta(4)*psi_simple_3;

%w_simple = theta(1)*psi_simple_0 + theta(2)*psi_simple_1 + theta(3)*psi_simple_2 + theta(4)*psi_simple_3;
w_simple = theta(1)*psi_simple_0 + theta(2)*psi_simple_1 + theta(3)*psi_simple_2 ;
%%
%plot of gaussians and resulting mexhat
filename = 'C:\Documents and Settings\lpolster\IDECorrData\src\matlab\scritps\Mexhat.pdf';
figure('units','centimeters','position',[0 5 plotwidth_fig2 plotheight_fig2],'filename',filename,...
   'papersize',[plotheight_fig2, plotwidth_fig2],'paperorientation','landscape','renderer','painters')  

plot(x,w_simple, 'black','LineWidth',3)
set(gca,'FontSize', 15)
xlabel('Space','FontSize', 15)
ylabel('Connectivity Strength','FontSize', 15)
hold on
plot(x,psi_simple_0_scaled, 'green')
hold on
plot(x,psi_simple_1_scaled,'red')
legend()
hold on
plot(x,psi_simple_2_scaled, 'blue')
% hold on
% plot(x,psi_simple_3_scaled,'magenta')
%legend('connectivity kernel','short-range excitation','mid-range inhibition','long-range excitation','anisotropic kernel')
legend('Mexican hat','short-range excitation','mid-range inhibition','long-range excitation')
% %%
% psi_0_large = Define2DGaussian(0,0, sigma_psi(1)^2, 0,2*NPoints-1,2*SpaceMin,2*SpaceMax);
% psi_1_large = Define2DGaussian(0,0, sigma_psi(2)^2, 0,2*NPoints-1,2*SpaceMin,2*SpaceMax);
% psi_2_large = Define2DGaussian(0,0, sigma_psi(3)^2, 0,2*NPoints-1,2*SpaceMin,2*SpaceMax);
% 
% w_large = theta(1)*psi_0_large + theta(2)*psi_1_large + theta(3)*psi_2_large;       % the large kernel %kernel k

filename = 'C:\Documents and Settings\lpolster\IDECorrData\src\matlab\scritps\KernelPlot.pdf'; % need to add
figure('units','centimeters','position',[10 5 plotwidth_fig3 plotheight_fig3],'filename',filename,...
   'papersize',[plotheight_fig3, plotwidth_fig3],'paperorientation','landscape','renderer','painters')  

cmax = 25;        % for plotting anisotropic: 10
cmin = -10.0;
imagesc(r,r,w,[cmin,cmax])
%title('Spatial Anisotropic Kernel','FontSize', 25)
xlabel('Space','FontSize', 15)
ylabel('Space','FontSize', 15)
xlim([-10,10])
ylim([-10,10])
set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10],'FontSize', 15)
axis square
axis xy
colorbar('FontSize', 15)

%% sensor parameters
% ~~~~~~~~~~~~~~~
sensor_index = 1:3:NPoints;                     % sets the location of the sensors
NSensors_xy = length(sensor_index);     % number of sensors in the x OR y plane
NSensors = NSensors_xy^2;                   % total number of sensors
sigma_m = 0.9;              % sensor width
m = Define2DGaussian(0, 0, sigma_m^2, 0, NPoints, SpaceMin, SpaceMax);% sensor kernel
m_large = Define2DGaussian(0, 0, sigma_m^2, 0, 2*NPoints-1,2*SpaceMin,2*SpaceMax);
delta_y = r(sensor_index(2)) - r(sensor_index(1));          % spacing between sensors in mm

%% observation noise characteristics
% ~~~~~~~~~~~~~~~~~~~~
sigma_varepsilon = 0.1;                                                                     
% Sigma_varepsilon = sigma_varepsilon^2*eye(NSensors);        % observation covariance matrix
Sigma_varepsilon = 0.1*eye(NSensors);        % observation covariance matrix

% create the observation noise
varepsilon = mvnrnd(zeros(1,NSensors),Sigma_varepsilon,T);
    
%% sigmoid parameters
% ~~~~~~~~~~~~~~~~ 
% f_max = 1;                  % maximum firing rate
varsigma = 0.56;         % sigmoid slope
v_0 = 1.8;                    % firing threshold, value from Wendling 2002??

%% synaptic kernel parameter
% ~~~~~~~~~~~~~~~~~~~~
tau = 0.01;                   % synaptic time constant
zeta = 1/tau;                 % inverse synaptic time constant
xi = 1-Ts*zeta;           % coefficient for the discrete time model



%% disturbance paramters
% ~~~~~~~~~~~~~
sigma_gamma = 1.3;              % parameter for covariance of disturbance
gamma_weight = 0.1;            % variance of disturbance

SphericalBoundary               % run the script that generates the covariance matrix

% gamma = gamma_weight*Define2DGaussian(0,0, sigma_gamma^2, 0,NPoints,SpaceMin,SpaceMax);
e = mvnrnd(zeros(1,NPoints^2),Sigma_gamma,T);

%% set initial condition
% ~~~~~~~~~~~~~~~~
v = zeros(T,NPoints,NPoints);
y = zeros(T,NSensors_xy,NSensors_xy);

max_field_init = gamma_weight;
min_field_init = -gamma_weight;
InitialCondition = min_field_init + (max_field_init - min_field_init)*rand(NPoints,NPoints);
v(1,:,:) = InitialCondition;

v_temp = padarray(InitialCondition,size(InitialCondition),'circular');
y_full = conv2(v_temp,m_large,'valid')*Delta_squared;      % the full field filtered by the larger sensor kernel
y_full = y_full(2:end-1,2:end-1); % trim to correct size

varepsilon_t = reshape(varepsilon(1,:,:), NSensors_xy, NSensors_xy);
y(1,:,:) = y_full(sensor_index,sensor_index) + varepsilon_t;

R_yy = zeros(T,2*NSensors_xy-1,2*NSensors_xy-1);  %initialize autocorrelation matrix
R_varepsilon_varepsilon = zeros(T,2*NSensors_xy-1,2*NSensors_xy-1);


% fighandle = figure;
%% Generate data
for t=1:T-1                                                                                                                                                                          
    
    y_t = squeeze(y(t,:,:));    
    v_t = squeeze(v(t,:,:));
    
    R_yy(t,:,:) = xcorr2(y_t);  %autocorrelation of observation function

    % calculate the firing rate     
    f = 1./(1+exp( varsigma*(v_0-v_t) ));           % calc firing rate using sigmoid
    f = padarray(f,size(f),'circular');
    g = conv2(f,w_large,'valid')*Delta_squared;   
    g = g(2:end-1,2:end-1);

    % conv firing rate with spatial kernel
    e_t = reshape(e(t,:,:),NPoints,NPoints);        % the disturbance
    v_t_plus1 = Ts*g + xi*v_t + e_t;  % update field %equation 3.10
    v(t+1,:,:) = v_t_plus1;

%     imagesc(v_t_plus1)
%     set(gca,'FontSize', 20)
%     colorbar ('FontSize', 20)
%     drawnow

    % filter field with sensor kernel and get observations
    v_t_plus1 = padarray(v_t_plus1,size(v_t_plus1),'circular');
    y_full = conv2(v_t_plus1,m_large,'valid') * Delta_squared;                  % field filtered by observation kernel at t+1
    y_full = y_full(2:end-1,2:end-1);

    varepsilon_tplus1 = reshape(varepsilon(t+1,:,:),NSensors_xy,NSensors_xy);
    y_tplus1 = y_full(sensor_index,sensor_index) + varepsilon_tplus1;                 % discretize sensor spacing equation 3.12
    y(t+1,:,:) = y_tplus1;
    
    R_varepsilon_varepsilon(t,:,:) = xcorr2(varepsilon_tplus1,varepsilon_tplus1);
    
   
    
    %% Plots of generated data
%     if ~(exist('Plots','dir') == 7) % folder name
%         mkdir('Plots');
%     end
%     
%     filenamebase = 'Plots/generated_Data_'; 
%     
%     if t == 500 || t == 501 ||  t == 502 || t == 503 ||  t == 504
%         
%         save_idx = 0;
%         tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.fig');
%         while exist(strcat(tempfilename),'file') == 2;
%             tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.fig');
%             save_idx = save_idx +1;
%         end
%         saveas(fighandle, tempfilename, 'fig');
%         
%         save_idx = 0;
%         tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.jpg');
%         while exist(strcat(tempfilename),'file') == 2;
%             tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.jpg');
%             save_idx = save_idx +1;
%         end
%         print(fighandle, '-djpeg', tempfilename);
%         
%         save_idx = 0;
%         tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.eps');
%         while exist(strcat(tempfilename),'file') == 2;
%             tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.eps');
%             save_idx = save_idx +1;
%         end
%         print(fighandle, '-depsc', tempfilename);
%         
%         save_idx = 0;
%         tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.pdf');
%         while exist(strcat(tempfilename),'file') == 2;
%             tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.pdf');
%             save_idx = save_idx +1;
%         end
%         print(fighandle, '-dpdf', tempfilename);
% 
%         clear tempfilename save_idx
       
     % change the size of the plot or resolution: http://www.mathworks.com.au/help/matlab/ref/print.html
     
        
 end

%%
%draw function at different points of time

% if t == 500 
% a=v_t_plus1(1:40,1:40);
% end
% if t == 501
% b=v_t_plus1(1:40,1:40);
% end
% if t == 502 
% c=v_t_plus1(1:40,1:40);
% end
% if t == 503 
% d=v_t_plus1(1:40,1:40);
% end

% end
% figure; hold on;
% 
% ha=pcolor(a);
% hb=pcolor(b);
% hc=pcolor(c);
% hd=pcolor(d);
% % he=pcolor(e);
% % hf=pcolor(f);
% shading flat
% shading interp
% 
% 
% set(hb,'zdata',0*b+2)
% set(hc,'zdata',0*c+4)
% set(hd,'zdata',0*d+6)
% axis([1 40 1 40])
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[]) 
% set(gca,'ZTickLabel',[])
% xlabel('Space','FontSize', 20)
% ylabel('Space','FontSize', 20)
% zlabel('Time','FontSize', 20)
% title('Evolution of the Spatial Field','FontSize', 25)
%%


%  mean_Noise = squeeze(mean(R_varepsilon_varepsilon(500:end-1,:,:),1));
%  R_corr = squeeze(mean(R_yy(500:end-1,:,:),1));
%  F_R_corr = fft2(R_corr - mean_Noise);
%     
% %     plot(min(F_R_corr))
% %     plot ((sigma_varepsilon)^2)
%     
%     diff = min((abs(F_R_corr)))-((sigma_varepsilon)^2);
%     figure
%     plot(diff)
%     figure
%     hist(diff)
toc

%%
% tempxxx = zeros(T,41,41);
% for n=1:T
%     tempxxx(n,:,:) = cov(squeeze(v(n,:,:)));
% 
% end
% figure,imagesc(squeeze(mean(tempxxx,1)))
% figure,imagesc(Sigma_gamma)
