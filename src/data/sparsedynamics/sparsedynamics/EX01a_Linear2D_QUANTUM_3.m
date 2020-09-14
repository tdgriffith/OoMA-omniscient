% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

clear all, close all, clc
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
% figpath = '../figures/';
% addpath('./utils');

%% generate Data
n=2;
basis=gen_SA_basis(n);
A = zeros(n,n);
for i=1:length(basis)
    A=A+(basis{i}*-rand());
end
%A=[-1, 2-3j;2+3j -4]
A=[-1,3j;-3j,-2];
A=A*-j; %make it quantum baby
rhs = @(x)A*x;   % ODE right hand side
dt=.01;
tspan=[0:dt:10];   % time span
x0=rand(n,1)+j*rand(n,1); %initial conditions
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
[t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate

%% compute Derivative 
eps = 0.0;      % noise strength
% for i=1:length(x)
%     dx(i,:) = A*x(i,:)';
% end

% dx_real=-1*real(gradient(x',dt)');
% dx_imag=imag(gradient(x',dt)');
% dx_split=dx_real+(dx_imag*j);
% dx=dx_split
dx=gradient(x',dt)';
%dx = dx  + eps*randn(size(dx));   % add noise
dx2=reshape(dx.',1,[])'
%% Custom Library
Theta_cust=[];
for i = 1:length(basis)
    x_sigma = x*(basis{i}*j);
    x_sigma = reshape(x_sigma.',1,[])';
    Theta_cust = [Theta_cust,x_sigma];
end
%Theta_cust=Theta_cust*-j;

%% compute Sparse regression: sequential least squares
Xi = [real(Theta_cust);imag(Theta_cust)]\[real(dx2);imag(dx2)]; %Really not sure here. Something messed up with imag
%Xi(4)=-Xi(4); %wtf, think its ordered wrong???
%Xi = imag(Theta_cust)\imag(dx2);
A_solve=zeros(n,n);
for i = 1:length(basis)
    A_solve=A_solve+(Xi(i)*basis{i}*-j);
end
A_solve
    

% lambda = 0.05;      % lambda is our sparsification knob.
% Xi = sparsifyDynamics(Theta,dx,lambda,n)

%% integrate true and identified systems
[tA,xA]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
rhs2 = @(x)A_solve*x;   % ODE right hand side
[tB,xB]=ode45(@(t,x)rhs2(x),tspan,x0,options);  % integrate

%% FIGURES!!
A
A_solve
error=norm(A-A_solve)

figure
subplot(1,2,1)
plot(tA,xA,'LineWidth',1.5)
hold on
plot(tB,xB,'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
title('Real. Comp. A')

subplot(1,2,2)

plot(tA,imag(xA),'LineWidth',1.5)
hold on
plot(tB,imag(xB),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
title('Imag. Comp. A')
sgtitle(['Simulation Results True A vs. SysID A: Error ',num2str(error)])

