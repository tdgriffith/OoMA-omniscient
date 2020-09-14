clear all, close all
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
%% Import Data Restructuring
channel_locs=readtable('deap_imag_conv.csv');
%% Grab Data: Using the just H Model here
n=2;
basis=gen_SA_basis(n);

subject = 4;
if subject <= 9
    load_name1=['s0',num2str(subject),'.mat']
else
    load_name1=['s',num2str(subject),'.mat']
end
subject1=subject;

s01=load(load_name1);
disp(subject);

trial = 18;
Y1=s01.data(trial,1:32,:);
Y1=squeeze(Y1);
Y1=Y1';
real_channs=channel_locs.Real_DEAP_(1:n);
imag_channs=channel_locs.IMAG_DEAP_(1:n);
Y1_real=Y1(:,real_channs);
Y1_imag=Y1(:,imag_channs)*1j;
Y1=Y1_real+Y1_imag;
extension='.png';
fs=128;
dt=1/fs;
tspan=linspace(0,63,8064);
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

%n_offdiag=(((n*n)-n)/2)-1;
%Xi(end-n_offdiag:end)=Xi(end-n_offdiag:end)*-1;
A_solve=zeros(n,n);
for i = 1:length(basis)
    if ~isreal(basis{i})
        Xi(i)=-Xi(i)
    end
    A_solve=A_solve+(Xi(i)*basis{i});
end
A_solve

%error=norm(A-A_solve) how to error data data

%[tA,xA]=ode45(@(t,x)rhs(x),tspan,x0,options);   % true model
rhs2 = @(x)A_solve*x;   % ODE right hand side

options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
x0=Y1(1,:);
[tB,xB]=ode45(@(t,x)rhs2(x),tspan,x0,options);  % integrate
recon_error=norm(Y1-xB,'Fro');
figure
subplot(1,2,1)
plot(tspan(1:128),Y1(1:128,:),'LineWidth',1.5)
hold on
plot(tB(1:128),xB(1:128,:),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
title('Real. Comp. A')

subplot(1,2,2)

plot(tspan(1:128),imag(Y1(1:128,:)),'LineWidth',1.5)
hold on
plot(tspan(1:128),imag(xB(1:128,:)),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
title('Imag. Comp. A')
sgtitle(['Real Data SA Test, Fro Norm: ', num2str(recon_error)])