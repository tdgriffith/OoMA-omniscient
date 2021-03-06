%% OMA with Complexity Plots for Average Over Subjects
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
onlineratings=gen_onlineratings();
%% Setup
Phi_cov_out=[];
Phi_data_out=[];
Phi_NeXT_out=[];
Phi_n4sid_out=[];
fn_cov_out=[];
fn_data_out=[];
fn_NeXT_out=[];
fn_n4sid_out=[];
A_cov_out=[];
A_data_out=[];
A_NeXT_out=[];
A_n4_out=[];
C_cov_out=[];
C_data_out=[];
C_NeXT_out=[];
C_n4_out=[];
extension='.png';


for subject = 2:2
    if subject <= 9
        load_name1=['S0',num2str(subject),'_restingPre_EC.mat']
    else
        load_name1=['S',num2str(subject),'_restingPre_EC.mat']
    end
    subject1=subject;
    s01=load(load_name1);
    
    % Loop Over Trials
    
    % Trial Data
    Y1=s01.dataRest(1:30,:);
    Y1=bandpass(detrend(Y1'),[3 45],fs);
    t = linspace(0,60*2.5,length(Y1));
    fs=2048;
    T=1/fs;
    
    % OMA Covariance Algorithm
    order = 60;
    s = 2*order;

    [A_cov,C_cov,G_cov,R0_cov,S_cov] = ssicov(Y1,order,s);
    opt_order=sum(S_cov/max(S_cov)>.001);
    
    eig(A_cov{opt_order})
    err = [0.01,0.05,0.98];
    [IDs_cov] = plotstab(A_cov,C_cov,Y1,T,[],err);
    title(['Stab. Diagram for Resting Subject ', num2str(subject), ', Trial ' num2str(trial)])
    filename=['export/resting/S' num2str(subject),'_rest', extension]
    saveas(gcf,filename)
    
    [fn_cov,zeta_cov,Phi_cov] = modalparams(A_cov,C_cov,T);
    norm_cov=abs(max(Phi_cov{opt_order}, [], 'all'));
    Phi_cov{opt_order}=Phi_cov{opt_order}/norm_cov;
    fn_cov_out{subject}=fn_cov{opt_order};
    Phi_cov_out{subject}=Phi_cov{opt_order};
    A_cov_out{subject}=A_cov{opt_order};
    C_cov_out{subject}=C_cov{opt_order};
    
    
    %n4sid
    data1 = iddata(Y1,[],T);
    sys1 = n4sid(data1,opt_order);
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
    [IDs_cov] = plotstab(A_data,C_data,Y1,T,[],err);
    
    
    [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T);
    Phi_data=Phi_data';
    norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
    Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest
    Phi_data_out{subject}=Phi_data{opt_order};
    fn_data_out{subject}=fn_data{opt_order};
    A_data_out{subject}=A_data{opt_order};
    C_data_out{subject}=C_data{opt_order};
    
    %         %NeXT-ERA
    [NeXT]=NExTFERA(Y1',5,2000,4,0.1,fs,800,200,opt_order,10,1);
    
    
    [fn_NeXT,zeta_NeXT,Phi_NeXT] = modalparams(NeXT.Matrices.A,NeXT.Matrices.C,T);
    Phi_NeXT=Phi_NeXT{1};
    norm_NeXT=abs(max(Phi_NeXT, [], 'all'));
    Phi_NeXT=Phi_NeXT/norm_NeXT;
    Phi_NeXT_out{subject}=Phi_NeXT;
    fn_NeXT_out{subject}=fn_NeXT{1};
    A_NeXT_out{subject}=NeXT.Matrices.A;
    C_NeXT_out{subject}=NeXT.Matrices.C;
    
    
    
%     figure
%     set(gcf,'units','points','position',[500,-300,700,700])
%     
%     subplot(2,2,1)
%     
%     for i=1:size(Phi_cov{opt_order},2)
%         hidden_arrow = compass(1,0);
%         hidden_arrow.Color = 'none';
%         hold on
%         compass(Phi_cov{opt_order}(:,i),colors(i))
%         hold on
%     end
%     title('OMA-Covar')
%     xlim([-1.25 1.25])
%     ylim([-1.25 1.25])
%     grid on
%     
%     
%     % OMA-data plot
%     subplot(2,2,2)
%     for i=1:size(Phi_data{opt_order},2)
%         hidden_arrow = compass(1,0);
%         hidden_arrow.Color = 'none';
%         hold on
%         compass(Phi_data{opt_order}(:,i),colors(i))
%         hold on
%     end
%     title('OMA-data')
%     xlim([-1.25 1.25])
%     ylim([-1.25 1.25])
%     grid on
%     
%     subplot(2,2,3)
%     for i=1:size(Phi_n4,2)
%         hidden_arrow = compass(1,0);
%         hidden_arrow.Color = 'none';
%         hold on
%         compass(Phi_n4(:,i),colors(i))
%         hold on
%     end
%     title('n4sid')
%     xlim([-1.25 1.25])
%     ylim([-1.25 1.25])
%     grid on
%     
%     subplot(2,2,4)
%     for i=1:size(Phi_NeXT,2)
%         hidden_arrow = compass(1,0);
%         hidden_arrow.Color = 'none';
%         hold on
%         compass(Phi_NeXT(:,i),colors(i))
%         hold on
%     end
%     title('NeXT')
%     xlim([-1.25 1.25])
%     ylim([-1.25 1.25])
%     grid on
%     
%     sgtitle(join((['Eigenvector Complexity Plots for Resting Subject: ' num2str(subject)])))
%     filename=['export/resting/S' num2str(subject),'rest' ,extension]
%     
%     saveas(gcf,filename)
%     close all
    
    figure
    set(gcf,'units','points','position',[500,-200,700,500])
    zeta_map=transpose(zeta_data{opt_order}/max(zeta_data{opt_order}));
    fn_map=transpose(fn_data{opt_order}/max(fn_data{opt_order}));
    h_indmap=heatmap([real(Phi_data{opt_order});imag(Phi_data{opt_order});zeta_map;fn_map]);
    h_indmap.Colormap=parula;
    h_indmap.ColorbarVisible=0;
    
    %sgtitle(join((['Heatmaps (Averages of A) for Emotion: ' onlineratings(trial)])))
    filename2=['export/heatmaps/spis_resting/S' num2str(subject),'T' num2str(trial), extension]
    saveas(gcf,filename2)
    %close all
end
