% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

close all
clear all
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
load_mat='s08.mat';
s01=load(load_mat);
trial=20;
fs=128;
dt=1/fs;
%% generate Data
polyorder = 5;  % search space up to fifth order polynomials
usesine = 0;    % no trig functions
n = 2;          % 2D system
sigma1=[0 1;1 0];
sigma2 = [0 -j;j 0];
sigma3 = [1 0;0 -1];
sigma4=eye(2);

x =s01.data(trial,[10,20],:);
x=squeeze(x);
tspan=linspace(0,size(x,2)/fs,size(x,2));
dx = gradient(x) ./ gradient(t);
%% compute Derivative 
eps = .0005;      % noise strength
dx = dx;  %+ eps*randn(size(dx));   % add noise
dx2=reshape(dx.',1,[])';
x=x'
%% Custom Library
x_sigma4 = x*sigma4;
x_sigma4 = real(reshape(x_sigma4.',1,[])');
x_sigma1 = x*sigma1;
x_sigma1 = real(reshape(x_sigma1.',1,[])');
% x_sigma2 = x*[0 -j;j 0];
% x_sigma2 = reshape(x_sigma2.',1,[])';
x_sigma2 = x*[0 -1;1 0];
x_sigma2 = imag(reshape(x_sigma2.',1,[])');
x_sigma3 = x*sigma3;
x_sigma3 = real(reshape(x_sigma3.',1,[])');
Theta_cust=[x_sigma1,x_sigma2,x_sigma3,x_sigma4];
%Theta_cust=[x_sigma1,x_sigma2,x_sigma4];

%% pool Data  (i.e., build library of nonlinear time series)

m = size(Theta_cust,2);

%% compute Sparse regression: sequential least squares
Xi = [real(Theta_cust);imag(Theta_cust)]\[real(dx2);imag(dx2)];
A_solve = (Xi(1)*sigma1)+(Xi(2)*sigma2)+(Xi(3)*sigma3)+(Xi(4)*sigma4);
%A_solve = (Xi(1)*sigma1)+(Xi(2)*sigma2)+(Xi(3)*sigma4);

% lambda = 0.05;      % lambda is our sparsification knob.
% Xi = sparsifyDynamics(Theta,dx,lambda,n)

%% integrate true and identified systems
[tA,xA]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
rhs2 = @(x)A_solve*x;   % ODE right hand side
[tB,xB]=ode45(@(t,x)rhs2(x),t,x0,options);  % integrate

%% FIGURES!!
figure
dtA = [0; diff(tA)];
plot(xA(:,1),xA(:,2),'r','LineWidth',1.5);
hold on
dtB = [0; diff(tB)];
plot(xB(:,1),xB(:,2),'k--','LineWidth',1.2);
xlabel('x_1','FontSize',13)
ylabel('x_2','FontSize',13)
l1 = legend('True','Identified');

set(gcf,'units','points','position',[10,10,750,400])
subplot(1,2,1)
plot(tA,xA(:,1),'r','LineWidth',1.5)
hold on
plot(tA,xA(:,2),'b-','LineWidth',1.5)
plot(tB(1:10:end),xB(1:10:end,1),'k--','LineWidth',1.2)
hold on
plot(tB(1:10:end),xB(1:10:end,2),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, x_k')
legend('True x_1','True x_2','Identified')
title('Real Comp. A')

subplot(1,2,2)

plot(tA,imag(xA(:,1)),'r','LineWidth',1.5)
hold on
plot(tA,imag(xA(:,2)),'b-','LineWidth',1.5)
plot(tB(1:10:end),imag(xB(1:10:end,1)),'k--','LineWidth',1.2)
hold on
plot(tB(1:10:end),imag(xB(1:10:end,2)),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, x_k')
legend('True x_1','True x_2','Identified')
title('Imag. Comp. A')
sgtitle('Simulation Results True A vs. SysID A')

A
A_solve