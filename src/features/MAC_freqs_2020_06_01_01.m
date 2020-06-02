%% OMA with Complexity Plots for Long Time Series from Original Data set H5 FORMAT
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
onlineratings=gen_onlineratings();
%% Setup
Phi_cov_out{1,32}=[];
Phi_data_out{1,32}=[];
Phi_NeXT_out{1,32}=[];
Phi_n4sid_out{1,32}=[];
fn_cov_out{1,32}=[];
fn_data_out{1,32}=[];
fn_NeXT_out{1,32}=[];
fn_n4sid_out{1,32}=[];
A_cov_out{1,32}=[];
A_data_out{1,32}=[];
A_NeXT_out{1,32}=[];
A_n4_out{1,32}=[];
C_cov_out{1,32}=[];
C_data_out{1,32}=[];
C_NeXT_out{1,32}=[];
C_n4_out{1,32}=[];
fs1=128;
fs2=120;
T1=1/fs1;
T2=1/fs2;

parfor subject = 1:32
    if subject <= 9
        load_name1=['s0',num2str(subject),'_',num2str(fs1),'.h5']
        load_name2=['s0',num2str(subject),'_',num2str(fs2),'.h5']
    else
        load_name1=['s',num2str(subject),'_',num2str(fs1),'.h5']
        load_name2=['s',num2str(subject),'_',num2str(fs2),'.h5']
    end
    subject1=subject;
    s01=h5read(load_name1,['/df_',num2str(fs1),'/block0_values']);
    s01_detrend=s01';
    s01_detrend=detrend(s01_detrend);
    s01=s01_detrend';
    
    s02=h5read(load_name2,['/df_',num2str(fs2),'/block0_values']);
    s02_detrend=s02';
    s02_detrend=detrend(s02_detrend);
    s02=s02_detrend';

    basename=['S' num2str(subject1)];
    extension='.png';
    
    %mkdir(['T',num2str(trial)])
    
    % Loop Over Trials
    
    % Trial Data
    channels=[3,4,31,27,28];
    Y1=s01(channels, length(s01)*.15:length(s01)*.3);
    Y2=s02(channels,length(s02)*.15:length(s02)*.3);
%     Y1=reshape(permute(Y1, [3 1 2]), [], size(Y1,2));
    Y1=Y1';
    Y2=Y2';
%     Y_train=Y1(:,1:322560*.8);
%     Y_test=Y1(:,322560*.8:322560);


    
    % OMA Covariance Algorithm
    order = 40;
    s = 2*order;
    opt_order=30;
    [A_cov,C_cov,G_cov,R0_cov,S_cov] = ssicov(Y1,order,s);
    [A_cov2,C_cov2,G_cov2,R0_cov2,S_cov2] = ssicov(Y2,order,s);
    figure
    bar(S_cov)
    title('Singular Values for OMA-COV')
    filename=['export/',num2str(fs1),'hz_filter/S' num2str(subject),'sigma',extension]
    %saveas(gcf,filename)
    close all
    
    %eig(A_cov{opt_order})
    err = [0.01,0.05,0.98];
%     [IDs_cov] = plotstab(A_cov,C_cov,Y1,T1,[],err);
%     [IDs_cov2] = plotstab(A_cov2,C_cov2,Y2,T2,[],err);
% 
%     filename=['export/',num2str(fs1),'hz_filter/S' num2str(subject),'stabCOV',extension]
    %saveas(gcf,filename)
    %close all
    
    [fn_cov,zeta_cov,Phi_cov] = modalparams(A_cov,C_cov,T1);
    [fn_cov2,zeta_cov2,Phi_cov2] = modalparams(A_cov2,C_cov2,T2);
    norm_cov=abs(max(Phi_cov{opt_order}, [], 'all'));
    norm_cov2=abs(max(Phi_cov2{opt_order}, [], 'all'));
    Phi_cov{opt_order}=Phi_cov{opt_order}/norm_cov;
    Phi_cov2{opt_order}=Phi_cov2{opt_order}/norm_cov2;
%     fn_cov_out{subject}=fn_cov{opt_order};
%     fn_cov_out2{subject}=fn_cov2{opt_order};
%     Phi_cov_out{subject}=Phi_cov{opt_order};
%     A_cov_out{subject}=A_cov{opt_order};
%     C_cov_out{subject}=C_cov{opt_order};
    
    
    
    %n4sid: Running really slow with 512hz
%     data1 = iddata(Y1,[],T);
%     %         data2 = iddata(Y2',[],T);
%     sys1 = n4sid(data1,opt_order);
%     %         sys2 = n4sid(data2,opt_order);
%     [A_n4,~,C_n4,~] = ssdata(sys1);
%     
%     
%     [fn_n4,zeta_n4,Phi_n4] = modalparams(A_n4,C_n4,T);
%     Phi_n4=Phi_n4{1};
%     norm_n4=abs(max(Phi_n4, [], 'all'));
%     Phi_n4=Phi_n4/norm_n4;
%     Phi_n4_out{subject}=Phi_n4;
%     fn_n4sid_out{subject}=fn_n4{1};
%     A_n4_out{subject}=A_n4;
%     C_n4_out{subject}=C_n4;
    
    % OMA-data
    [A_data,C_data,G_data,R0_data] = ssidata(Y1,order,s);
    [A_data2,C_data2,G_data2,R0_data2] = ssidata(Y2,order,s);
    
    
    [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T1);
    [fn_data2,zeta_data2,Phi_data2] = modalparams(A_data2,C_data2,T2);
%     [IDs_cov] = plotstab(A_data,C_data,Y1,T1,[],err);
%     filename=['export/',num2str(fs1),'hz_filter/S' num2str(subject),'stabDATA',extension];
    %saveas(gcf,filename)
    %close all
    
    Phi_data=Phi_data';
    norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
    Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest
    Phi_data2=Phi_data2';
    norm_data2=abs(max(Phi_data2{opt_order}, [], 'all')); %determine max value
    Phi_data2{opt_order}=Phi_data2{opt_order}/norm_data2; %normalize the one of interest
    
    %MAC plots
    figure
    mm_cov=macmatrix(Phi_cov{10},Phi_cov2{10});
    macplot(mm_cov)
    filename=['export/MAC_freqs/',num2str(fs1),'hz_MAC/S' num2str(subject),'_OMAcov',extension];
    saveas(gcf,filename)
    close all
    
    figure
    mm_data=macmatrix(Phi_data{10},Phi_data2{10});
    macplot(mm_data)
    filename=['export/MAC_freqs/',num2str(fs1),'hz_MAC/S' num2str(subject),'_OMAdata',extension];
    saveas(gcf,filename)
    close all


end
% Save numbers at end
    
    







    