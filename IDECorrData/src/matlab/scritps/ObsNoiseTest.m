clc
clear
close all

Ts = 1e-3;          % sampling period (s)
T = 50000;            % maximum time (ms)

%% spatial parameters
% ~~~~~~~~~~~
Delta = 0.5;                          % space step for the spatial discretisation
Delta_squared = Delta^2;
SpaceMax = 10;                    % maximum space in mm
SpaceMin = -SpaceMax;         % minimum space in mm
NPoints = (SpaceMax-SpaceMin)/Delta+1;
NPoints_total = NPoints^2;
r = linspace(SpaceMin,SpaceMax,NPoints);      % define space

sensor_index = 1:3:41;
NSensors_xy = length(sensor_index);
NSensors = NSensors_xy^2;
sigma_m = 0.9;                                                         % sensor width
m = Define2DGaussian(0, 0, sigma_m^2, 0, NPoints, SpaceMin, SpaceMax);

%% observation noise characteristics
% ~~~~~~~~~~~~~~~~~~~~
sigma_varepsilon = 0.01; %0.1                                                                    
Sigma_varepsilon = sigma_varepsilon^2*eye(NSensors);        % observation covariance matrix

% create the observation noise
varepsilon = mvnrnd(zeros(1,NSensors),Sigma_varepsilon,T);



R_noise = zeros(T,2*NSensors_xy-1,2*NSensors_xy-1);

for t=1:T-1
    Noise = reshape(varepsilon(t+1,:,:),NSensors_xy,NSensors_xy);
    R_noise(t+1,:,:) = xcorr2(Noise);
end

sigma_varepsilon^2 * NSensors

R_Noise_mean = squeeze(mean(R_noise,1));
R_Noise_mean(14,14)

figure,imagesc(R_Noise_mean),colorbar
