%% OMA with Complexity Plots for Average Over Subjects
set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
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

for trial= 1:1
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

        % Loop Over Trials

        % Trial Data
        Y1=s01.data(trial,[3,4,31,27,28],:);
        Y1=squeeze(Y1);
        %Y2=s01.data(trial+1,[3,4,31,27,28],:);
        %Y2=squeeze(Y2);
        t = linspace(0,63,length(Y1));
        fs=128;
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
        %norm_n4=abs(max(Phi_n4, [], 'all'));
        %Phi_n4=Phi_n4/norm_n4;
        Phi_n4_out{subject}=Phi_n4;
        fn_n4sid_out{subject}=fn_n4{1};
        A_n4_out{subject}=A_n4;
        C_n4_out{subject}=C_n4;

        % OMA-data
        [A_data,C_data,G_data,R0_data] = ssidata(Y1,order,s);


        [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T);
        Phi_data=Phi_data';
        %norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
        %Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest
        Phi_data_out{subject}=Phi_data{opt_order};
        fn_data_out{subject}=fn_data{opt_order};
        A_data_out{subject}=A_data{opt_order};
        C_data_out{subject}=C_data{opt_order};

        %         %NeXT-ERA
        [NeXT]=NExTFERA(Y1,5,2000,4,0.1,fs,800,200,opt_order,10,1);


        [fn_NeXT,zeta_NeXT,Phi_NeXT] = modalparams(NeXT.Matrices.A,NeXT.Matrices.C,T);
        Phi_NeXT=Phi_NeXT{1};
        %norm_NeXT=abs(max(Phi_NeXT, [], 'all'));
        %Phi_NeXT=Phi_NeXT/norm_NeXT;
        Phi_NeXT_out{subject}=Phi_NeXT;
        fn_NeXT_out{subject}=fn_NeXT{1}; 
        A_NeXT_out{subject}=NeXT.Matrices.A;
        C_NeXT_out{subject}=NeXT.Matrices.C;
    end
    % ReShape bois
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
    % Averages
    fn_cov_3d=cat(3,fn_cov_out{:});
    fn_cov_mean=mean(fn_cov_3d,3);
    fn_cov_covar=cov(squeeze(fn_cov_3d));

    Phi_cov_3d=cat(3,Phi_cov_out{:});
    Phi_cov_mean=mean(Phi_cov_3d,3);
    norm_cov=max(abs(Phi_cov_mean), [], 'all');
    Phi_cov_mean=Phi_cov_mean/norm_cov;

    Phi_data_3d=cat(3,Phi_data_out{:});
    Phi_data_mean=mean(Phi_data_3d,3);
    norm_data=max(abs(Phi_data_mean), [], 'all'); %determine max value
    Phi_data_mean=Phi_data_mean/norm_data; %normalize the one of interest

    Phi_n4_3d=cat(3,Phi_n4_out{:});
    Phi_n4_mean=mean(Phi_n4_3d,3);
    norm_n4=max(abs(Phi_n4_mean), [], 'all');
    Phi_n4_mean=Phi_n4_mean/norm_n4;

    Phi_NeXT_3d=cat(3,Phi_NeXT_out{:});
    Phi_NeXT_mean=mean(Phi_NeXT_3d,3);
    norm_NeXT=max(abs(Phi_NeXT_mean), [], 'all');
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
    
    filename=['export/averagePhi/T' num2str(trial),'avg' ,extension]
    
    saveas(gcf,filename)
    close all
    
    % Scatter Plots
    x=[1 1 1 1 1 1;2 2 2 2 2 2;3 3 3 3 3 3;4 4 4 4 4 4;5 5 5 5 5 5];
    [theta_cov,rho_cov]=cart2pol(real(Phi_cov_3d),imag(Phi_cov_3d));
    [theta_data,rho_data]=cart2pol(real(Phi_data_3d),imag(Phi_data_3d));
    [theta_n4,rho_n4]=cart2pol(real(Phi_n4_3d),imag(Phi_n4_3d));
    [theta_NeXT,rho_NeXT]=cart2pol(real(Phi_NeXT_3d),imag(Phi_NeXT_3d));
    
    theta_cov_mean=mean(theta_cov,3);
    rho_cov_mean=mean(rho_cov,3);
    rho_cov_norm=max(abs(rho_cov_mean), [], 'all');
    rho_cov_mean=rho_cov_mean./rho_cov_norm;
    N=32;
    SEM_theta_cov=std(theta_cov,0,3)/sqrt(N);
    SEM_rho_cov=std(rho_cov,0,3)/sqrt(N);
    
    
    theta_data_mean=mean(theta_data,3);
    rho_data_mean=mean(rho_data,3);
    rho_data_norm=max(abs(rho_data_mean), [], 'all');
    rho_data_mean=rho_data_mean/rho_data_norm;
    SEM_theta_data=std(theta_data,0,3)/sqrt(N);
    SEM_rho_data=std(rho_data,0,3)/sqrt(N);
    
    theta_n4_mean=mean(theta_n4,3);
    rho_n4_mean=mean(rho_n4,3);
    rho_n4_norm=max(abs(rho_n4_mean), [], 'all');
    rho_n4_mean=rho_n4_mean/rho_n4_norm;
    SEM_theta_n4=std(theta_n4,0,3)/sqrt(N);
    SEM_rho_n4=std(rho_n4,0,3)/sqrt(N);
    
    
    
    theta_NeXT_mean=mean(theta_NeXT,3);
    rho_NeXT_mean=mean(rho_NeXT,3);
    rho_NeXT_norm=max(abs(rho_NeXT_mean), [], 'all');
    rho_NeXT_mean=rho_NeXT_mean/rho_NeXT_norm;
    SEM_theta_NeXT=std(theta_cov,0,3)/sqrt(N);
    SEM_rho_NeXT=std(rho_cov,0,3)/sqrt(N);
    
    figure
    set(gcf,'units','points','position',[500,-300,1000,500])
    subplot(4,2,1)
    errorbar(x,theta_cov_mean,SEM_theta_cov,'-s')
    title('Angle Components for OMA-Cov')
    subplot(4,2,2)
    errorbar(x,rho_cov_mean,SEM_rho_cov,'-s')
    title('Radial Component for OMA-Cov')
    
    subplot(4,2,3)
    errorbar(x,theta_data_mean,SEM_theta_data,'-s')
    title('Angle Components for OMA-Data')
    subplot(4,2,4)
    errorbar(x,rho_data_mean,SEM_rho_data,'-s')
    title('Radial Component for OMA-Data')
    
    subplot(4,2,5)
    errorbar(x,theta_n4_mean,SEM_theta_n4,'-s')
    title('Angle Components for OMA-n4sid')
    subplot(4,2,6)
    errorbar(x,rho_n4_mean,SEM_rho_n4,'-s')
    title('Radial Component for OMA-n4sid')
    
    subplot(4,2,7)
    errorbar(x,theta_NeXT_mean,SEM_theta_NeXT,'-s')
    title('Angle Components for OMA-NeXT')
    subplot(4,2,8)
    errorbar(x,rho_NeXT_mean,SEM_rho_NeXT,'-s')
    title('Radial Component for OMA-NeXT')
    sgtitle('Standard Error of the Mean for Eigenvector Components')
    
    filename=['export/averagePhi/T' num2str(trial),'comps' ,extension]
    saveas(gcf,filename)
    close all
    
    %% Average of A's and Herms
    % n4
    A_n4_out_3d=cat(3,A_n4_out{:});
    A_n4_average=mean(A_n4_out_3d,3);
    A_n4_herm=A_n4_average+A_n4_average';
    C_n4_out_3d=cat(3,C_n4_out{:});
    C_n4_average=mean(C_n4_out_3d,3);
    [fn_n4_2, zeta_n4_2, Phi_n4_2]=modalparams(A_n4_average,C_n4_average,T);
    [fn_n4_3, zeta_n4_3, Phi_n4_3]=modalparams(A_n4_herm,C_n4_average,T);
    Phi_n4_2=Phi_n4_2{1};
    norm_n4_2=max(abs(Phi_n4_2), [],'all');
    Phi_n4_2=Phi_n4_2/norm_n4_2;
    Phi_n4_3=Phi_n4_3{1};
    norm_n4_3=max(abs(Phi_n4_3), [],'all');
    Phi_n4_3=Phi_n4_3/norm_n4_3;
    
    % NeXT
    A_NeXT_out_3d=cat(3,A_NeXT_out{:});
    A_NeXT_average=mean(A_NeXT_out_3d,3);
    A_NeXT_herm=A_NeXT_average+A_NeXT_average';
    C_NeXT_out_3d=cat(3,C_NeXT_out{:});
    C_NeXT_average=mean(C_NeXT_out_3d,3);
    [fn_NeXT_2, zeta_NeXT_2, Phi_NeXT_2]=modalparams(A_NeXT_average,C_NeXT_average,T);
    [fn_NeXT_3, zeta_NeXT_3, Phi_NeXT_3]=modalparams(A_NeXT_herm,C_NeXT_average,T);
    Phi_NeXT_2=Phi_NeXT_2{1};
    norm_NeXT_2=max(abs(Phi_NeXT_2), [],'all');
    Phi_NeXT_2=Phi_NeXT_2/norm_NeXT_2;
    Phi_NeXT_3=Phi_NeXT_3{1};
    norm_NeXT_3=max(abs(Phi_NeXT_3), [],'all');
    Phi_NeXT_3=Phi_NeXT_3/norm_NeXT_3;
    
    % OMA-COV
    A_cov_out_3d=cat(3,A_cov_out{:});
    A_cov_average=mean(A_cov_out_3d,3);
    A_cov_herm=A_cov_average+A_cov_average';
    C_cov_out_3d=cat(3,C_cov_out{:});
    C_cov_average=mean(C_cov_out_3d,3);
    [fn_cov_2, zeta_cov_2, Phi_cov_2]=modalparams(A_cov_average,C_cov_average,T);
    [fn_cov_3, zeta_cov_3, Phi_cov_3]=modalparams(A_cov_herm,C_cov_average,T);
    Phi_cov_2=Phi_cov_2{1};
    norm_cov_2=max(abs(Phi_cov_2), [],'all');
    Phi_cov_2=Phi_cov_2/norm_cov_2;
    Phi_cov_3=Phi_cov_3{1};
    norm_cov_3=max(abs(Phi_cov_3), [],'all');
    Phi_cov_3=Phi_cov_3/norm_cov_3;
    
    % data
    A_data_out_3d=cat(3,A_data_out{:});
    A_data_average=mean(A_data_out_3d,3);
    A_data_herm=A_data_average+A_data_average';
    C_data_out_3d=cat(3,C_data_out{:});
    C_data_average=mean(C_data_out_3d,3);
    [fn_data_2, zeta_data_2, Phi_data_2]=modalparams(A_data_average,C_data_average,T);
    [fn_data_3, zeta_data_3, Phi_data_3]=modalparams(A_data_herm,C_data_average,T);
    Phi_data_2=Phi_data_2{1};
    norm_data_2=max(abs(Phi_data_2), [],'all');
    Phi_data_2=Phi_data_2/norm_data_2;
    Phi_data_3=Phi_data_3{1};
    norm_data_3=max(abs(Phi_data_3), [],'all');
    Phi_data_3=Phi_data_3/norm_data_3;
    
    figure
    set(gcf,'units','points','position',[500,-300,700,700])
    
    subplot(2,2,1)
    
    for i=1:length(Phi_cov_2)
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_cov_2(:,i),colors(i))
        hold on
    end
    title('OMA-Covar')
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])
    grid on
    
    
    % OMA-data plot
    subplot(2,2,2)
    
    for i=1:length(Phi_data_2)
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_data_2(:,i),colors(i))
        hold on
    end
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])
    grid on
    title('OMA-Data')
    
    % Plot n4sid
    subplot(2,2,3)
    
    for i=1:length(Phi_n4_2(1,:))
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_n4_2(:,i),colors(i))
        hold on
    end
    title('n4sid')
    
    %Plot NeXT-ERA
    subplot(2,2,4)
    
    
    for i=1:length(Phi_NeXT_2(1,:))
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_NeXT_2(:,i),colors(i))
        hold on
    end
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])
    grid on
    title('NeXT-ERA')
    
    
    sgtitle(join((['Eigenvector Complexity Plots (Averages of A) for Emotion: ' onlineratings(trial)])))
    
    filename=['export/averageA/T' num2str(trial),'avgofA' ,extension]
    
    saveas(gcf,filename)
    close all
    
    figure
    set(gcf,'units','points','position',[500,-300,700,700])
    
    subplot(2,2,1)
    
    for i=1:length(Phi_cov_3)
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_cov_3(:,i),colors(i))
        hold on
    end
    title('OMA-Covar')
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])
    grid on
    
    
    % OMA-data plot
    subplot(2,2,2)
    
    for i=1:length(Phi_data_3)
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_data_3(:,i),colors(i))
        hold on
    end
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])
    grid on
    title('OMA-Data')
    
    % Plot n4sid
    subplot(2,2,3)
    
    for i=1:length(Phi_n4_3(1,:))
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_n4_3(:,i),colors(i))
        hold on
    end
    title('n4sid')
    
    %Plot NeXT-ERA
    subplot(2,2,4)
    
    
    for i=1:length(Phi_NeXT_3(1,:))
        hidden_arrow = compass(1,0);
        hidden_arrow.Color = 'none';
        hold on
        compass(Phi_NeXT_3(:,i),colors(i))
        hold on
    end
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])
    grid on
    title('NeXT-ERA')
    
    
    sgtitle(join((['Eigenvector Complexity Plots (Averages of A) for Emotion: ' onlineratings(trial)])))
    
    filename=['export/hermA/T' num2str(trial),'hermofA' ,extension]
    
    saveas(gcf,filename)
    close all
    
    
    

end





    