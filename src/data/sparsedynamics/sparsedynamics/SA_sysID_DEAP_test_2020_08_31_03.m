clear all, close all
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
%% Import Data Restructuring
channel_locs=readtable('deap_imag_conv.csv');
%% Grab Data: Using the just H Model here
n=2;
basis=gen_iH_basis(n);

subject = 1;
if subject <= 9
    load_name1=['s0',num2str(subject),'.mat']
else
    load_name1=['s',num2str(subject),'.mat']
end
subject1=subject;

s01=load(load_name1);
disp(subject);

trial = 12;
Y1=s01.data(trial,1:32,:);
Y1=squeeze(Y1);
Y1=Y1';

extension='.png';
fs=128;
dt=1/fs;
tspan=linspace(0,63,8064);

order = 2*n;
s = 2*order;
opt_order=n;
[A_data,C_data,G_data,R0_data] = ssidata(Y1,order,s);
Ah=fwht(A_data{opt_order});
[fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,dt);

A=Ah*-1j;
rhs = @(x)A*x;   % ODE right hand side
tspan=[0:.01:10];   % time span
x0=rand(n,1)+1j*rand(n,1); %initial conditions
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
[t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate
Y1=x;

%% compute Derivative 
eps = 0.0;      % noise strength
dx=gradient(Y1.',dt).';

%% Custom Library
Theta_cust=[];
for i = 1:length(basis)
    x_sigma = Y1*(basis{i});
    Theta_cust = [Theta_cust,x_sigma];
end
Theta_cust=reshape(Theta_cust,[],length(basis));
Theta=[real(Theta_cust);imag(Theta_cust)];

%% Estimate State Matrix
dx3=[real(dx(:));imag(dx(:))];
Xi=mldivide(Theta,dx3);

lambda=max(abs(Xi))*0.2;
smallinds = (abs(Xi)<lambda);   % find small coefficients
biginds = ~smallinds;
Theta=Theta(:,biginds);
Xi=mldivide(Theta,dx3);
basis(smallinds) = [];

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
        Xi(i)=-Xi(i);
    end
    A_solve=A_solve+(Xi(i)*basis{i});
end
A_solve;
A
A_solve=gen_quant_A(A)
%error=norm(A-A_solve) how to error data data

%[tA,xA]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
rhs2 = @(x)A_solve*x;   % ODE right hand side

options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
x0=Y1(1,:);
[tB,xB]=ode45(@(t,x)rhs2(x),tspan,x0,options);  % integrate
recon_error=norm(Y1-xB,'Fro');
figure
subplot(1,2,1)
plot(tspan(1:512),Y1(1:512,:),'LineWidth',1.5)
hold on
plot(tB(1:512),xB(1:512,:),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
title('Real. Comp. A')

subplot(1,2,2)

plot(tspan(1:512),imag(Y1(1:512,:)),'LineWidth',1.5)
hold on
plot(tspan(1:512),imag(xB(1:512,:)),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
title('Imag. Comp. A')
sgtitle(['Real Data SA Test, Fro Norm: ', num2str(recon_error)])