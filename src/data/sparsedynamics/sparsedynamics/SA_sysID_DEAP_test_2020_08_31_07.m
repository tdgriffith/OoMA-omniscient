clear all, close all
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
%% Import Data Restructuring
%channel_locs=readtable('deap_imag_conv.csv');
%% Grab Data: Trying DCT INSTEAD OF OMA 
max_chans=16;
subject = 1;
if subject <= 9
    load_name1=['s0',num2str(subject),'.mat']
else
    load_name1=['s',num2str(subject),'.mat']
end
subject1=subject;

s01=load(load_name1);
disp(subject);

trial = 4;
Y1=s01.data(trial,1:max_chans,:);
Y1=squeeze(Y1);
Y1=Y1';
B_ts=[];
1*4096/64;
for i = 1:128
    B=fft2(Y1(i:i+7,1:8)); %or dct2
    B_vec=B(:);
    B_vec=B_vec./norm(B_vec); %RISKY BOIS HERE
    B_ts=[B_ts,B_vec];
end
B_ts=B_ts';
figure
subplot(1,2,1)
plot(real(B_ts(1:100,1:10)))
subplot(1,2,2)
plot(imag(B_ts(1:100,1:10)))

extension='.png';
fs=128;
dt=1/fs;

%% compute Derivative 
eps = 0.0;      % noise strength
dx=gradient(B_ts.',dt).';

basis=gen_iH_basis(size(B_ts,2));

Theta_cust=[];
for i = 1:length(basis)
    x_sigma = (basis{i})*x.';
    Theta_cust = [Theta_cust,x_sigma.'];
    disp(i/length(basis))
end

Theta_cust=reshape(Theta_cust,[],length(basis));
Theta=[real(Theta_cust);imag(Theta_cust)];

%% Estimate State Matrix
dx3=[real(dx(:));imag(dx(:))];
Xi=mldivide(Theta,dx3);

lambda=5;
smallinds = (abs(Xi)<lambda);   % find small coefficients
biginds = ~smallinds;
Theta=Theta(:,biginds);
Xi=mldivide(Theta,dx3);
basis(smallinds) = [];


A_solve=zeros(size(B_ts,2),size(B_ts,2));
for i = 1:length(basis)
    if isreal(basis{i})
        Xi(i)=-Xi(i);
    end
    A_solve=A_solve+(Xi(i)*basis{i});
    disp(i/length(basis))
end
A_solve;

rhs2 = @(x)A_solve*x;   % ODE right hand side
tspan=[0:dt:10];

ch=8;
ch2=8;
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,size(B_ts,2)));

[tB,xB]=ode45(@(t,x)rhs2(x),tspan,B_ts(1,:),options);  % integrate

figure
subplot(1,2,1)
plot(tspan(1:200),B_ts(1:200,ch:ch2),'LineWidth',1.5)
hold on
plot(tB(1:200),xB(1:200,ch:ch2),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
ylim([-1,1])
title('Real. Comp. A')

subplot(1,2,2)

plot(tspan(1:200),imag(B_ts(1:200,ch:ch2)),'LineWidth',1.5)
hold on
plot(tB(1:200),imag(xB(1:200,ch:ch2)),'k--','LineWidth',1.2)
xlabel('Time')
ylabel('State, $x_k$')
ylim([-1,1])
title('Imag. Comp. A')
recon_error=norm(xB(1:200)-B_ts(1:200));
sgtitle(['Real Data SA Test, Fro Norm: ', num2str(recon_error)])