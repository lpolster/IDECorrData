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
FS_Label = 8;          % point / fontsize for the axis label
FS_Tick = 8;                % fontsize for the ticks
MS = 10;                     % marker size
LW = 1;
plotwidth = 4.4;        % cm
plotheight = 3.5;

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
T = 1000;              % maximum time (ms)          % 2 seconds = 100 seconds

%% spatial kernel parameters
% ~~~~~~~~~~~~~~

theta(1) = 100;%80.0;           % local kernel amplitude
theta(2) = -80;             % surround kernel amplitude
theta(3) = 5;               % lateral kernel amplitude

sigma_psi(1) = 1.8;     % local kernel width
sigma_psi(2) = 2.4;     % surround kernel width
sigma_psi(3) = 6;       % lateral kernel width

psi_0 = Define2DGaussian(0,0, sigma_psi(1)^2, 0,NPoints,SpaceMin,SpaceMax);
psi_1 = Define2DGaussian(0,0, sigma_psi(2)^2, 0,NPoints,SpaceMin,SpaceMax);
psi_2 = Define2DGaussian(0,0, sigma_psi(3)^2, 0,NPoints,SpaceMin,SpaceMax);

psi_0_scaled = theta(1)*psi_0;
psi_1_scaled = theta(2)*psi_1;
psi_2_scaled = theta(3)*psi_2;

w = theta(1)*psi_0 + theta(2)*psi_1 + theta(3)*psi_2;       % the kernel

% plotting 3 scaled gaussians and resulting mexican hat
% filename = 'C:\Documents and Settings\lpolster\IDECorrData\src\matlab\scritps\4KernelPlot.pdf';
% figure('filename',filename)  
% subplot(221)
% surf(w)
% title('Spatial Kernel')
% colorbar
% subplot(222)
% surf(psi_0_scaled)
% title('Local Kernel')
% colorbar
% subplot(223)
% surf(psi_1_scaled)
% title('Surround Kernel')
% colorbar
% subplot(224)
% surf(psi_2_scaled)
% title('Lateral Kernel')
% colorbar

psi_0_large = Define2DGaussian(0,0, sigma_psi(1)^2, 0,2*NPoints-1,2*SpaceMin,2*SpaceMax);
psi_1_large = Define2DGaussian(0,0, sigma_psi(2)^2, 0,2*NPoints-1,2*SpaceMin,2*SpaceMax);
psi_2_large = Define2DGaussian(0,0, sigma_psi(3)^2, 0,2*NPoints-1,2*SpaceMin,2*SpaceMax);

w_large = theta(1)*psi_0_large + theta(2)*psi_1_large + theta(3)*psi_2_large;       % the large kernel

filename = 'C:\Documents and Settings\lpolster\IDECorrData\src\matlab\scritps\KernelPlot.pdf'; % need to add
%figure('units','centimeters','position',[0 0 plotwidth plotheight],'filename',filename,...
%    'papersize',[plotheight, plotwidth],'paperorientation','landscape','renderer','painters')  

% cmax = 25.5;        % for plotting
% cmin = -10.0;
% imagesc(r,r,w,[cmin,cmax])
% title('Spatial Kernel')
% xlabel('Space','fontsize',FS_Label)
% ylabel('Space','fontsize',FS_Label)
% xlim([-10,10])
% ylim([-10,10])
% set(gca,'xtick',[-10 0 10],'ytick',[-10 0 10],'fontsize',FS_Tick)
% axis square
% axis xy
% colorbar

%% sensor parameters
% ~~~~~~~~~~~~~~~
sensor_index = 1:3:NPoints;                     % sets the location of the sensors
NSensors_xy = length(sensor_index);     % number of sensors in the x OR y plane
NSensors = NSensors_xy^2;                   % total number of sensors
sigma_m = 0.9;              % sensor width
m = Define2DGaussian(0, 0, sigma_m^2, 0, NPoints, SpaceMin, SpaceMax);      % sensor kernel
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

fighandle = figure;
%% Generate data
for t=1:T-1                                                                                                                                                                          
    
    y_t = squeeze(y(t,:,:));    
    v_t = squeeze(v(t,:,:));

    % calculate the firing rate     
    f = 1./(1+exp( varsigma*(v_0-v_t) ));           % calc firing rate using sigmoid

    f = padarray(f,size(f),'circular');
    g = conv2(f,w_large,'valid')*Delta_squared;   
    g = g(2:end-1,2:end-1);

    % conv firing rate with spatial kernel
    e_t = reshape(e(t,:,:),NPoints,NPoints);        % the disturbance
    v_t_plus1 = Ts*g + xi*v_t + e_t;  % update field 
    v(t+1,:,:) = v_t_plus1;
   
imagesc(v_t_plus1)
colorbar
drawnow

    % filter field with sensor kernel and get observations
    v_t_plus1 = padarray(v_t_plus1,size(v_t_plus1),'circular');
    y_full = conv2(v_t_plus1,m_large,'valid') * Delta_squared;                  % field filtered by observation kernel at t+1
    y_full = y_full(2:end-1,2:end-1);

    varepsilon_tplus1 = reshape(varepsilon(t+1,:,:),NSensors_xy,NSensors_xy);
    y_tplus1 = y_full(sensor_index,sensor_index) + varepsilon_tplus1;                 % discretize sensor spacing
    y(t+1,:,:) = y_tplus1;
    
    
    
    if ~(exist('Plots','dir') == 7) % folder name
        mkdir('Plots');
    end
    
    filenamebase = 'Plots/generated_Data_'; 
    
    if t == 500 || t == 501 ||  t == 502 || t == 503 ||  t == 504
        
        save_idx = 0;
        tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.fig');
        while exist(strcat(tempfilename),'file') == 2;
            tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.fig');
            save_idx = save_idx +1;
        end
        saveas(fighandle, tempfilename, 'fig');
        
        save_idx = 0;
        tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.jpg');
        while exist(strcat(tempfilename),'file') == 2;
            tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.jpg');
            save_idx = save_idx +1;
        end
        print(fighandle, '-djpeg', tempfilename);
        
        save_idx = 0;
        tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.eps');
        while exist(strcat(tempfilename),'file') == 2;
            tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.eps');
            save_idx = save_idx +1;
        end
        print(fighandle, '-depsc', tempfilename);
        
        save_idx = 0;
        tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.pdf');
        while exist(strcat(tempfilename),'file') == 2;
            tempfilename = strcat(filenamebase,num2str(t),'_Run',num2str(save_idx),'.pdf');
            save_idx = save_idx +1;
        end
        print(fighandle, '-dpdf', tempfilename);

        clear tempfilename save_idx
       
     % change the size of the plot or resolution: http://www.mathworks.com.au/help/matlab/ref/print.html
     
        
    end
    
end

toc

%%
% tempxxx = zeros(T,41,41);
% for n=1:T
%     tempxxx(n,:,:) = cov(squeeze(v(n,:,:)));
% 
% end
% figure,imagesc(squeeze(mean(tempxxx,1)))
% figure,imagesc(Sigma_gamma)