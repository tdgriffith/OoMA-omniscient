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
fs=512;
T=1/fs;

parfor subject = 7:7
    if subject <= 9
        load_name1=['s0',num2str(subject),'_',num2str(fs),'.h5']
    else
        load_name1=['s',num2str(subject),'_',num2str(fs),'.h5']
    end
    subject1=subject;
    s01=h5read(load_name1,['/df_',num2str(fs),'/block0_values']);
    s01_detrend=s01';
    s01_detrend=detrend(s01_detrend);
    s01=s01_detrend';

    basename=['S' num2str(subject1)];
    extension='.png';
    
    %mkdir(['T',num2str(trial)])
    
    % Loop Over Trials
    
    % Trial Data
    channels=[3,4,31,27,28];
    Y1=s01(:,:);
%     Y1=reshape(permute(Y1, [3 1 2]), [], size(Y1,2));
    Y1=Y1';
%     Y_train=Y1(:,1:322560*.8);
%     Y_test=Y1(:,322560*.8:322560);


    
    % OMA Covariance Algorithm
    order = 40;
    s = 2*order;
    opt_order=30;
    [A_cov,C_cov,G_cov,R0_cov,S_cov] = ssicov(Y1,order,s);
    figure
    bar(S_cov)
    title('Singular Values for OMA-COV')
    filename=['export/',num2str(fs),'hz_filter/S' num2str(subject),'sigma',extension]
    %saveas(gcf,filename)
    %close all
    
    %eig(A_cov{opt_order})
    err = [0.01,0.05,0.98];
    [IDs_cov] = plotstab(A_cov,C_cov,Y1,T,[],err);
    filename=['export/',num2str(fs),'hz_filter/S' num2str(subject),'stabCOV',extension]
    %saveas(gcf,filename)
    %close all
    
    [fn_cov,zeta_cov,Phi_cov] = modalparams(A_cov,C_cov,T);
    norm_cov=abs(max(Phi_cov{opt_order}, [], 'all'));
    Phi_cov{opt_order}=Phi_cov{opt_order}/norm_cov;
    fn_cov_out{subject}=fn_cov{opt_order};
    Phi_cov_out{subject}=Phi_cov{opt_order};
    A_cov_out{subject}=A_cov{opt_order};
    C_cov_out{subject}=C_cov{opt_order};
    
    
    
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
    
    
    [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T);
    [IDs_cov] = plotstab(A_data,C_data,Y1,T,[],err);
    filename=['export/',num2str(fs),'hz_filter/S' num2str(subject),'stabDATA',extension];
    %saveas(gcf,filename)
    %close all
    
    Phi_data=Phi_data';
    norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
    Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest
    Phi_data_out{subject}=Phi_data{opt_order};
    fn_data_out{subject}=fn_data{opt_order};
    A_data_out{subject}=A_data{opt_order};
    C_data_out{subject}=C_data{opt_order};
    
    %         %NeXT-ERA
    n_windows=(length(Y1)/20000)-8;
    [NeXT]=NExTFERA(Y1',1,20000,n_windows,0.1,fs,800,200,2*opt_order,10,1);
    
    
    [fn_NeXT,zeta_NeXT,Phi_NeXT] = modalparams(NeXT.Matrices.A,NeXT.Matrices.C,T);
    Phi_NeXT=Phi_NeXT{1};
    norm_NeXT=abs(max(Phi_NeXT, [], 'all'));
    Phi_NeXT=Phi_NeXT/norm_NeXT;
    Phi_NeXT_out{subject}=Phi_NeXT;
    fn_NeXT_out{subject}=fn_NeXT{1};
    A_NeXT_out{subject}=NeXT.Matrices.A;
    C_NeXT_out{subject}=NeXT.Matrices.C;
    
    figure
    set(gcf,'units','points','position',[500,-300,700,700])
    
    subplot(2,2,1)
    
    for i=1:length(Phi_cov{opt_order})
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_cov{opt_order}(:,i),colors(i))
        hold on
    end
    title('OMA-Covar')
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])
    grid on
%     
%     
%     % OMA-data plot
    subplot(2,2,2)
    
    for i=1:length(Phi_data{opt_order}) 
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_data{opt_order}(:,i),colors(i))
        hold on
    end
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])
    grid on
    title('OMA-Data')
%     
%     % Plot n4sid
    subplot(2,2,3)
    
    for i=1:1 %length(Phi_n4) FIX ME DUMMY
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
%         hold on
%         compass(Phi_n4(:,i),colors(i))
%         hold on
    end
    title('n4sid')
%     
%     %Plot NeXT-ERA
    subplot(2,2,4)
    
    
    for i=1:length(Phi_NeXT(1,:))
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_NeXT(:,i),colors(i))
        hold on
    end
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])
    grid on
    title('NeXT-ERA')
%     
%     
    sgtitle(join((['Eigenvector Complexity Plots for ',num2str(fs),' hz, Subject ' num2str(subject)])))
    
    filename=['export/',num2str(fs),'hz_filter/S' num2str(subject),'compPlot',extension]
    
    %saveas(gcf,filename)
    %close all
    
end
% Save numbers at end
    
    







    