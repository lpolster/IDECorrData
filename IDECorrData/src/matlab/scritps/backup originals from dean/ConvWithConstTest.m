% convolution with a constant test
clc
close all
clear

varsigma = 0.56;
v_0 = 1.8;

Delta = 0.5;                          % space step for the spatial discretisation
Delta_squared = Delta^2;
SpaceMax = 10;                    % maximum space in mm
SpaceMin = -SpaceMax;         % minimum space in mm
NPoints = (SpaceMax-SpaceMin)/Delta+1;

sigma_m = 0.9;        % sensor width
sensor_index = 1:3:41;
NSensors_xy = length(sensor_index);

m = Define2DGaussian(0, 0, sigma_m^2, 0, NPoints, SpaceMin, SpaceMax);

theta(1) = 100.0;           % local kernel amplitude
theta(2) = -80;             % surround kernel amplitude
theta(3) = 5;               % lateral kernel amplitude

sigma_psi(1) = 1.8;     % local kernel width
sigma_psi(2) = 2.4;     % surround kernel width
sigma_psi(3) = 6;       % lateral kernel width

psi_0 = Define2DGaussian(0,0, sigma_psi(1)^2, 0,NPoints,SpaceMin,SpaceMax);
psi_1 = Define2DGaussian(0,0, sigma_psi(2)^2, 0,NPoints,SpaceMin,SpaceMax);
psi_2 = Define2DGaussian(0,0, sigma_psi(3)^2, 0,NPoints,SpaceMin,SpaceMax);
w = theta(1)*psi_0 + theta(2)*psi_1 + theta(3)*psi_2;       % the kernel


constant_field1 = (2-varsigma*v_0)*ones(2*NPoints,2*NPoints);

TestResult1 = conv2(w,constant_field1)*Delta_squared;       % convolve w with the constant field

figure,imagesc(TestResult1),colorbar

TestResult1(ceil(size(TestResult1,1)/2),ceil(size(TestResult1,2)/2))
c_1 =  sum(sum(w))*Delta_squared * (2-varsigma*v_0)

%%
constant_field2 = c_1*ones(2*NPoints,2*NPoints);

TestResult2_temp = conv2(m,constant_field2,'same')*Delta_squared;
TestResult2 = TestResult2_temp(sensor_index,sensor_index);

figure,imagesc(TestResult2),colorbar

TestResult2(ceil(size(TestResult2,1)/2),ceil(size(TestResult2,2)/2))
c_2 = sum(sum(m)) * Delta_squared * c_1

%%
y = randn(NSensors_xy,NSensors_xy);

constant_field3 = c_2*ones(2*NPoints,2*NPoints);

TestResult3 = xcorr2(constant_field3,y);
figure,imagesc(TestResult3),colorbar

TestResult3(ceil(size(TestResult3,1)/2),ceil(size(TestResult3,2)/2))
c_3 = sum(sum(y))*c_2

