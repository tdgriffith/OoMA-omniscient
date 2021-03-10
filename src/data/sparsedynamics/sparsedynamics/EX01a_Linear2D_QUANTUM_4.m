% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

clear all, close all
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
% figpath = '../figures/';
% addpath('./utils');

%% generate Data
n=2;
basis=gen_iH_basis(n);
A = zeros(n,n);
for i=1:length(basis)
    A=A+(basis{i}*rand());
end
%A=[-1, 0;0 0];
A=[0, j;-j 0];
A=A*j;

rhs = @(x)A*x;   % ODE right hand side
dt=.01;
tspan=[0:dt:10];   % time span
x0=rand(n,1)+j*rand(n,1); %initial conditions
x0=[1,0];
%x0=[-1+1j;1-1j]
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
[t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate


%% compute Derivative 
eps = 0.0;      % noise strength
dx=gradient(x.',dt).';
dx=dx+eps*randn(size(dx));

%% Custom Library
Theta_cust=[];
for i = 1:length(basis)
    x_sigma = x*(basis{i});
    Theta_cust = [Theta_cust,x_sigma];
end
Theta_cust=reshape(Theta_cust,[],length(basis));
Theta=[real(Theta_cust);imag(Theta_cust)];

%% Estimate State Matrix
dx3=[real(dx(:));imag(dx(:))];
Xi=mldivide(Theta,dx3);

lambda=0.1;
smallinds = (abs(Xi)<lambda);   % find small coefficients
biginds = ~smallinds;
Theta=Theta(:,biginds);
Xi=mldivide(Theta,dx3);
basis(smallinds) = [];

%n_offdiag=(((n*n)-n)/2)-1;
%Xi(end-n_offdiag:end)=Xi(end-n_offdiag:end)*-1;
A_solve=zeros(n,n);
for i = 1:length(basis)
    if isreal(basis{i})
        Xi(i)=-Xi(i)
    end
    A_solve=A_solve+(Xi(i)*basis{i});
end
A_solve
A
error=norm(A-A_solve)

[tA,xA]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
rhs2 = @(x)A_solve*x;   % ODE right hand side
[tB,xB]=ode45(@(t,x)rhs2(x),tspan,x0,options);  % integrate
normx=norm(xA);
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