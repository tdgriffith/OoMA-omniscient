%% OMA with Complexity Plots for Average Over Subjects using power
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
%% Setup
Phi_cov_out=[];
Phi_data_out=[];
Phi_NeXT_out=[];
Phi_n4sid_out=[];
fn_cov_out=[];
fn_data_out=[];
fn_NeXT_out=[];
fn_n4sid_out=[];
for subject = 1:32
    if subject <= 9
        load_name1=['s0',num2str(subject),'.mat']
    else
        load_name1=['s',num2str(subject),'.mat']
    end
    subject1=subject;
    s01=load(load_name1);
    load_int2=32-subject+1
    if load_int2 <=9
        load_name2=['s0',num2str(load_int2),'.mat']
    else
        load_name2=['s',num2str(load_int2),'.mat']
    end
    subject2=load_int2;
    basename=['S' num2str(subject1)];
    extension='.png';
    %s02=load(load_name2)
    
    %mkdir(['T',num2str(trial)])
    
    %% Loop Over Trials
    
    % Trial Data
    channels=[3,4,31,27,28];
    Y1=s01.data(:,channels,:);
    Y1=squeeze(Y1);
    Y1=reshape(Y1,[1,5,8064*40]);
    
    Y1=squeeze(Y1);
    Y1=Y1';
    
    Pxx=cell(1,length(channels));
    F=cell(1,length(channels));
    tvec=cell(1,length(channels));
    for iii=1:length(channels)
        [Pxx{iii},F{iii},tvec{iii}]=compute_eeg_power(Y1(:,iii),128,63*40);
    end
    ind_delta=find(F{1}<=4);
    ind_theta=find(F{1}>=4 & F{1}<=8);
    ind_alpha=find(F{1}>=8 & F{1}<=12);
    ind_beta=find(F{1}>=12 & F{1}<=30);
    ind_all=find(F{1}<=30);
    mean_delta=cell(1,length(channels));
    mean_alpha=cell(1,length(channels));
    mean_beta=cell(1,length(channels));
    mean_theta=cell(1,length(channels));
    Y_final=[];
    for iii=1:length(channels)
        mean_alpha{iii}=mean(Pxx{iii}(ind_alpha,:));
        mean_beta{iii}=mean(Pxx{iii}(ind_beta,:));
        mean_theta{iii}=mean(Pxx{iii}(ind_theta,:));
        mean_delta{iii}=mean(Pxx{iii}(ind_delta,:));
        Y_final=[Y_final, mean_delta{iii}', mean_theta{iii}',mean_alpha{iii}', mean_beta{iii}'];
    end
%     figure
%     plot(tvec{1},Y_final)
%     xlabel('Time (s)')
%     ylabel('Power (dB)')
%     title('Power Time Series for All Channels')
%     t = linspace(0,63,length(Y1));
    fs=128;
    fs=367/63;
    T=1/fs;
    
    % OMA Covariance Algorithm
    order = 25;
    s = 2*order;
    opt_order=12;
    [A_cov,C_cov,G_cov,R0_cov] = ssicov(Y1,order,s);
    
    eig(A_cov{opt_order})
    err = [0.01,0.05,0.98];
    %[IDs_cov] = plotstab(A_cov,C_cov,Y,T,[],err);
    [fn_cov,zeta_cov,Phi_cov] = modalparams(A_cov,C_cov,T);
    %norm_cov=abs(max(Phi_cov{opt_order}, [], 'all'));
    %Phi_cov{opt_order}=Phi_cov{opt_order}/norm_cov;
    fn_cov_out{subject}=fn_cov{opt_order};
    Phi_cov_out{subject}=Phi_cov{opt_order};
    
    
    %n4sid
    data1 = iddata(Y1',[],T);
    %         data2 = iddata(Y2',[],T);
    sys1 = n4sid(data1,opt_order);
    %         sys2 = n4sid(data2,opt_order);
    [A_n4_1,~,C_n4_1,~] = ssdata(sys1);
    
    
    [fn_n4,zeta_n4,Phi_n4] = modalparams(A_n4_1,C_n4_1,T);
    Phi_n4=Phi_n4{1};
    %norm_n4=abs(max(Phi_n4, [], 'all'));
    %Phi_n4=Phi_n4/norm_n4;
    Phi_n4_out{subject}=Phi_n4;
    fn_n4sid_out{subject}=fn_n4{1};
    
    % OMA-data
    [A_data,C_data,G_data,R0_data] = ssidata(Y1,order,s);
    
    
    [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T);
    Phi_data=Phi_data';
    %norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
    %Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest
    Phi_data_out{subject}=Phi_data{opt_order};
    fn_data_out{subject}=fn_data{opt_order};
    
    %         %NeXT-ERA
    [NeXT]=NExTFERA(Y1,5,2000,4,0.1,fs,800,200,opt_order,10,1);
    
    
    [fn_NeXT,zeta_NeXT,Phi_NeXT] = modalparams(NeXT.Matrices.A,NeXT.Matrices.C,T);
    Phi_NeXT=Phi_NeXT{1};
    %norm_NeXT=abs(max(Phi_NeXT, [], 'all'));
    %Phi_NeXT=Phi_NeXT/norm_NeXT;
    Phi_NeXT_out{subject}=Phi_NeXT;
    fn_NeXT_out{subject}=fn_NeXT{1};
end
%% ReShape bois
%reshape so all cells are same dim
for i=1:32
    fn_cov_out{i}=fn_cov_out{i}(1:min(cellfun('size',fn_cov_out,1)));
    Phi_cov_out{i}=Phi_cov_out{i}(:,1:min(cellfun('size',fn_cov_out,1))) ;
end
for i=1:32
    fn_data_out{i}=fn_data_out{i}(1:min(cellfun('size',fn_data_out,1)));
    Phi_data_out{i}=Phi_data_out{i}(:,1:min(cellfun('size',fn_data_out,1)));
end
for i=1:32
    fn_NeXT_out{i}=fn_NeXT_out{i}(1:min(cellfun('size',fn_NeXT_out,1)));
    Phi_NeXT_out{i}=Phi_NeXT_out{i}(:,1:min(cellfun('size',fn_NeXT_out,1)));
end
for i=1:32
    fn_n4sid_out{i}=fn_n4sid_out{i}(1:min(cellfun('size',fn_n4sid_out,1)));
    Phi_n4_out{i}=Phi_n4_out{i}(:,1:min(cellfun('size',fn_n4sid_out,1)));
end
%% Averages
fn_cov_3d=cat(3,fn_cov_out{:});
fn_cov_mean=mean(fn_cov_3d,3);
fn_cov_covar=cov(squeeze(fn_cov_3d));

Phi_cov_3d=cat(3,Phi_cov_out{:});
Phi_cov_mean=mean(Phi_cov_3d,3);
norm_cov=abs(max(Phi_cov_mean, [], 'all'));
Phi_cov_mean=Phi_cov_mean/norm_cov;

Phi_data_3d=cat(3,Phi_data_out{:});
Phi_data_mean=mean(Phi_data_3d,3);
norm_data=abs(max(Phi_data_mean, [], 'all')); %determine max value
Phi_data_mean=Phi_data_mean/norm_data; %normalize the one of interest

Phi_n4_3d=cat(3,Phi_n4_out{:});
Phi_n4_mean=mean(Phi_n4_3d,3);
norm_n4=abs(max(Phi_n4_mean, [], 'all'));
Phi_n4_mean=Phi_n4_mean/norm_n4;

Phi_NeXT_3d=cat(3,Phi_NeXT_out{:});
Phi_NeXT_mean=mean(Phi_NeXT_3d,3);
norm_NeXT=abs(max(Phi_NeXT_mean, [], 'all'));
Phi_NeXT_mean=Phi_NeXT_mean/norm_NeXT;

figure
set(gcf,'units','points','position',[500,-300,700,700])

subplot(2,2,1)

for i=1:length(Phi_cov_mean)
    hidden_arrow = compass(1,0);
    hidden_arrow.Color = 'none';
    hold on
    compass(Phi_cov_mean(:,i),colors(i))
    hold on
end
title('OMA-Covar')
xlim([-1.25 1.25])
ylim([-1.25 1.25])
grid on


% OMA-data plot
subplot(2,2,2)

for i=1:length(Phi_data_mean)
    hidden_arrow = compass(1,0);
    hidden_arrow.Color = 'none';
    hold on
    compass(Phi_data_mean(:,i),colors(i))
    hold on
end
xlim([-1.25 1.25])
ylim([-1.25 1.25])
grid on
title('OMA-Data')

% Plot n4sid
subplot(2,2,3)

for i=1:length(Phi_n4_mean(1,:))
    hidden_arrow = compass(1,0);
    hidden_arrow.Color = 'none';
    hold on
    compass(Phi_n4_mean(:,i),colors(i))
    hold on
end
title('n4sid')

%Plot NeXT-ERA
subplot(2,2,4)


for i=1:length(Phi_NeXT_mean(1,:))
    hidden_arrow = compass(1,0);
    hidden_arrow.Color = 'none';
    hold on
    compass(Phi_NeXT_mean(:,i),colors(i))
    hold on
end
xlim([-1.25 1.25])
ylim([-1.25 1.25])
grid on
title('NeXT-ERA')


sgtitle(join((['Eigenvector Complexity Plots (Averages) for Emotion: ' onlineratings(trial)])))

filename=['T' num2str(trial),'avg' ,extension]

saveas(gcf,filename)
close all

%% Scatter Plots
[theta_cov,rho_cov]=cart2pol(real(Phi_cov_3d),imag(Phi_cov_3d));
[theta_data,rho_data]=cart2pol(real(Phi_data_3d),imag(Phi_data_3d));
[theta_n4,rho_n4]=cart2pol(real(Phi_n4_3d),imag(Phi_n4_3d));
[theta_NeXT,rho_NeXT]=cart2pol(real(Phi_NeXT_3d),imag(Phi_NeXT_3d));

theta_cov_mean=mean(theta_cov,3);
rho_cov_mean=mean(rho_cov,3);
rho_cov_norm=abs(max(rho_cov_mean, [], 'all'));
rho_cov_mean=rho_cov_mean./rho_cov_norm;

theta_data_mean=mean(theta_data,3);
rho_data_mean=mean(rho_data,3);
rho_data_norm=abs(max(rho_data_mean, [], 'all'));
rho_data_mean=rho_data_mean/rho_data_norm;

theta_n4_mean=mean(theta_n4,3);
rho_n4_mean=mean(rho_n4,3);
rho_n4_norm=abs(max(rho_n4_mean, [], 'all'));
rho_n4_mean=rho_n4_mean/rho_n4_norm;

theta_NeXT_mean=mean(theta_NeXT,3);
rho_NeXT_mean=mean(rho_NeXT,3);
rho_NeXT_norm=abs(max(rho_NeXT_mean, [], 'all'));
rho_NeXT_mean=rho_NeXT_mean/rho_NeXT_norm;





