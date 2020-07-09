set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
onlineratings=gen_onlineratings();
filename_out=[];
%% Setup: Heatmaps for Robots NOT NORMALIZED OMA DATA
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

    
    %mkdir(basename)

%% Loop Over Trials
    for trial = 1:40
        % Trial Data
        Y1=s01.data(trial,[3,4,31,27,28],:);
        Y1=squeeze(Y1);
        %Y2=s01.data(trial+1,[3,4,31,27,28],:);
        %Y2=squeeze(Y2);
        t = linspace(0,63,length(Y1));
        fs=128;
        T=1/fs;

        % OMA Covariance Algorithm
        order = 60;
        s = 2*order;
        opt_order=22;
        %[A_cov,C_cov,G_cov,R0_cov] = ssicov(Y1,order,s);
%         [A_cov2,C_cov2,G_cov2,R0_cov2] = ssicov(Y2,order,s); 
%         eig(A_cov{opt_order})
        err = [0.01,0.05,0.98];
%         %[IDs_cov] = plotstab(A_cov,C_cov,Y,T,[],err);
%         [fn_cov,zeta_cov,Phi_cov] = modalparams(A_cov,C_cov,T);
%         norm_cov=abs(max(Phi_cov{opt_order}, [], 'all'));
%         Phi_cov{opt_order}=Phi_cov{opt_order}/norm_cov;
        
%         [fn_cov2,zeta_cov2,Phi_cov2] = modalparams(A_cov2,C_cov2,T);
%         real_part=real(eig(A_cov{opt_order}));
%         imag_part=imag(eig(A_cov{opt_order}));
%         real_part2=real(eig(A_cov2{opt_order}));
%         imag_part2=imag(eig(A_cov2{opt_order}));

        %n4sid
%         data1 = iddata(Y1',[],T);
% %         data2 = iddata(Y2',[],T);
%         sys1 = n4sid(data1,12);
% %         sys2 = n4sid(data2,12);
%         [A_n4_1,~,C_n4_1,~] = ssdata(sys1);
% %         [A_n4_2,~,C_n4_2,~] = ssdata(sys2);
% %         real_part_n4=real(eig(A_n4_1));
% %         imag_part_n4=imag(eig(A_n4_1));
% %         real_part_n4_2=real(eig(A_n4_2));
% %         imag_part_n4_2=imag(eig(A_n4_2));
%         
%         [fn_n4,zeta_n4,Phi_n4] = modalparams(A_n4_1,C_n4_1,T);
%         Phi_n4=Phi_n4{1};
%         norm_n4=abs(max(Phi_n4, [], 'all'));
%         Phi_n4=Phi_n4/norm_n4;

        % OMA-data
        [A_data,C_data,G_data,R0_data] = ssidata(Y1,order,s);
        [IDs_data] = plotstab(A_data,C_data,Y1,T,[],err);
        title(['Stab. Diagram for Subject ', num2str(subject), ', Trial ' num2str(trial)])
        mkdir(['export/stab/order60/S' num2str(subject)])
        filename=['export/stab/order60/S' num2str(subject) '/S' num2str(subject),'T' num2str(trial), extension]
        saveas(gcf,filename)
        %close all

        
%         [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T);
%         Phi_data=Phi_data';
%         norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
%         %Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest
% 
% %         %NeXT-ERA
%         [NeXT]=NExTFERA(Y1,5,2000,4,0.1,fs,800,200,12,10,1);
% %         real_part_NeXT=real(eig(NeXT.Matrices.A));
% %         imag_part_NeXT=imag(eig(NeXT.Matrices.A));
%         
%         [fn_NeXT,zeta_NeXT,Phi_NeXT] = modalparams(NeXT.Matrices.A,NeXT.Matrices.C,T);
%         Phi_NeXT=Phi_NeXT{1};
%         norm_NeXT=abs(max(Phi_NeXT, [], 'all'));
%         Phi_NeXT=Phi_NeXT/norm_NeXT;
        

        
        
    end
end