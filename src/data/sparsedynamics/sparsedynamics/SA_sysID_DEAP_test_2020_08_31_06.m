clear all, close all
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
%% Import Data Restructuring
channel_locs=readtable('deap_imag_conv.csv');
%% Grab Data: Try Method from 113.080401
%Need to know C matrix so this is out dfadfasf
A=[-4 0;0 -4]*1j;
H=hadamard(2)./norm(hadamard(2));


rhs = @(x)A*x;   % ODE right hand side
tspan=[0:.01:6];   % time span
x0=H*[1;-1j]; %initial conditions
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,2));
[tA,xA]=ode45(@(t,x)rhs(x),tspan,x0,options);  % integrate
[A_data,C_data,G_data,R0_data] = ssidata(xA,2,4);
[b,a]=ss2tf(A_data{2},[],C_data{2},[])





%% compute Derivative 


rhs2 = @(x)A*x;   % ODE right hand side
x0=[1;-1j];
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,2));
[tB,xB]=ode45(@(t,x)rhs2(x),tspan,x0,options);  % integrate
xB=xB.';
for i = 1:length(xB)
    %xB(:,i)=ifwht(xB(:,i));
    xB(:,i)=H*xB(:,i);
end
xB=xB.';
recon_error=norm(xA-xB,'Fro');
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
sgtitle(['Real Data SA Test, Fro Norm: ', num2str(recon_error)])