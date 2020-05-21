 %% OMA with Complexity Plots for Long Time Series
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
onlineratings=gen_onlineratings();
%% Setup
% Phi_cov_out=[];
% Phi_data_out=[];
% Phi_NeXT_out=[];
% Phi_n4sid_out=[];
% fn_cov_out=[];
% fn_data_out=[];
% fn_NeXT_out=[];
% fn_n4sid_out=[];
% A_cov_out=[];
% A_data_out=[];
% A_NeXT_out=[];
% A_n4_out=[];
% C_cov_out=[];
% C_data_out=[];
% C_NeXT_out=[];
% C_n4_out=[];


for subject = 1:1
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
    
    % Loop Over Trials
    
    % Trial Data
    Y1=s01.data(:,[3,4,31,27,28],:);
    Y1=reshape(permute(Y1, [3 1 2]), [], size(Y1,2));
    Y1=Y1';
    Y_train=Y1(:,1:322560*.8);
    Y_test=Y1(:,322560*.8:322560);
    %Y2=s01.data(trial+1,[3,4,31,27,28],:);
    %Y2=squeeze(Y2);
    fs=128;
    T=1/fs;
    
    % OMA Covariance Algorithm
    order = 30;
    s = 2*order;
    opt_order=22;
    [A_cov,C_cov,G_cov,R0_cov] = ssicov(Y1,order,s);
    
    eig(A_cov{opt_order})
    err = [0.01,0.05,0.98];
    %[IDs_cov] = plotstab(A_cov,C_cov,Y,T,[],err);
    [fn_cov,zeta_cov,Phi_cov] = modalparams(A_cov,C_cov,T);
    norm_cov=abs(max(Phi_cov{opt_order}, [], 'all'));
    Phi_cov{opt_order}=Phi_cov{opt_order}/norm_cov;
    fn_cov_out{subject}=fn_cov{opt_order};
    Phi_cov_out{subject}=Phi_cov{opt_order};
    A_cov_out{subject}=A_cov{opt_order};
    C_cov_out{subject}=C_cov{opt_order};
    
    
    %n4sid
    data1 = iddata(Y1',[],T);
    %         data2 = iddata(Y2',[],T);
    sys1 = n4sid(data1,opt_order);
    %         sys2 = n4sid(data2,opt_order);
    [A_n4,~,C_n4,~] = ssdata(sys1);
    
    
    [fn_n4,zeta_n4,Phi_n4] = modalparams(A_n4,C_n4,T);
    Phi_n4=Phi_n4{1};
    norm_n4=abs(max(Phi_n4, [], 'all'));
    Phi_n4=Phi_n4/norm_n4;
    Phi_n4_out{subject}=Phi_n4;
    fn_n4sid_out{subject}=fn_n4{1};
    A_n4_out{subject}=A_n4;
    C_n4_out{subject}=C_n4;
    
    % OMA-data
    [A_data,C_data,G_data,R0_data] = ssidata(Y1,order,s);
    
    
    [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T);
    Phi_data=Phi_data';
    norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
    Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest
    Phi_data_out{subject}=Phi_data{opt_order};
    fn_data_out{subject}=fn_data{opt_order};
    A_data_out{subject}=A_data{opt_order};
    C_data_out{subject}=C_data{opt_order};
    
    %         %NeXT-ERA
    [NeXT]=NExTFERA(Y1,5,2000,4,0.1,fs,800,200,opt_order,10,1);
    
    
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
    
    for i=1:length(Phi_n4)
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_n4(:,i),colors(i))
        hold on
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
    sgtitle(join((['Eigenvector Complexity Plots for ALL Trials, Subject ' num2str(subject)])))
    
    filename=['export/longTS/S' num2str(subject),'avg' ,extension]
    
    saveas(gcf,filename)
    close all
end
% ReShape bois
%reshape so all cells are same dim
    
    







    





    