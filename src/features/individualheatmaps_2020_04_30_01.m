set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
onlineratings=gen_onlineratings();
filename_out=[];
%% Setup: Heatmaps for Robots NORMALIZED OMA DATA
for subject = 1:32 %parfor
    if subject <= 9
        load_name1=['s0',num2str(subject),'.mat']
    else
        load_name1=['s',num2str(subject),'.mat']
    end
    subject1=subject;
    s01=load(load_name1);
    basename=['S' num2str(subject1)];
    extension='.png';

    
    %mkdir(basename)

%% Loop Over Trials
    parfor trial = 1:40
        % Trial Data
        %Y1=s01.data(trial,[3,4,31,27,28],:);
        Y1=s01.data(trial,1:32,:);
        Y1=squeeze(Y1);
        t = linspace(0,63,length(Y1));
        fs=128;
        T=1/fs;
        extension='.png';

        % OMA Covariance Algorithm
        order = 40;
        s = 2*order;
        opt_order=40;
        [A_cov,C_cov,G_cov,R0_cov,S_cov] = ssicov(Y1,order,s);
%         [A_cov2,C_cov2,G_cov2,R0_cov2] = ssicov(Y2,order,s); 
        eig(A_cov{opt_order})
        err = [0.01,0.05,0.98];
        %[IDs_cov] = plotstab(A_cov,C_cov,Y1,T,[],err);
        [fn_cov,zeta_cov,Phi_cov] = modalparams(A_cov,C_cov,T);
        norm_cov=abs(max(Phi_cov{opt_order}, [], 'all'));
        Phi_cov{opt_order}=Phi_cov{opt_order}/norm_cov;
        
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
        err = [0.01,0.05,0.98];
        %[IDs_cov] = plotstab(A_data,C_data,Y1,T,[],err);
%         real_part_data=real(eig(A_data{opt_order}));
%         imag_part_data=imag(eig(A_data{opt_order}));
%         [A_data2,C_data2,G_data2,R0_data2] = ssidata(Y2,order,s);
%         real_part_data2=real(eig(A_data2{opt_order}));
%         imag_part_data2=imag(eig(A_data2{opt_order}));
        
        [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T);
        Phi_data=Phi_data';
        norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
        Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest

%         %NeXT-ERA
%         [NeXT]=NExTFERA(Y1,5,2000,4,0.1,fs,800,200,12,10,1);
% %         real_part_NeXT=real(eig(NeXT.Matrices.A));
% %         imag_part_NeXT=imag(eig(NeXT.Matrices.A));
%         
%         [fn_NeXT,zeta_NeXT,Phi_NeXT] = modalparams(NeXT.Matrices.A,NeXT.Matrices.C,T);
%         Phi_NeXT=Phi_NeXT{1};
%         norm_NeXT=abs(max(Phi_NeXT, [], 'all'));
%         Phi_NeXT=Phi_NeXT/norm_NeXT;
        
        
%         [NeXT2]=NExTFERA(Y2,5,2000,4,0.1,fs,800,200,12,10,1);
%         real_part_NeXT2=real(eig(NeXT2.Matrices.A));
%         imag_part_NeXT2=imag(eig(NeXT2.Matrices.A));

        


        
        %heatmaps for robuts
        figure
        set(gcf,'units','points','position',[500,-200,700,500])
        zeta_map=transpose(zeta_data{opt_order}/max(zeta_data{opt_order}));
        fn_map=transpose(fn_data{opt_order}/max(fn_data{opt_order}));
        h_indmap=heatmap([real(Phi_data{opt_order});imag(Phi_data{opt_order});zeta_map;fn_map]);
        h_indmap.Colormap=parula;
        h_indmap.ColorbarVisible=0;
        
        %sgtitle(join((['Heatmaps (Averages of A) for Emotion: ' onlineratings(trial)])))
        filename2=['export/heatmaps/robots_allch/S' num2str(subject),'T' num2str(trial), extension]
        saveas(gcf,filename2)
        close all
        
%         %heatmaps for people
%         figure
%         set(gcf,'units','points','position',[500,-200,700,500])
%         h_indmap=heatmap([real(Phi_data{opt_order});imag(Phi_data{opt_order});zeta_map;fn_map]);
%         h_indmap.Colormap=parula;
%         h_indmap.ColorbarVisible=0;
%         %h_indmap.XDisplayLabels={'Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Mode 5', 'Mode 6', 'Mode 7', 'Mode 8', 'Mode 9', 'Mode 10', 'Mode 11', 'Mode 12', 'Mode 13', 'Mode 14', 'Mode 15', 'Mode 16', 'Mode 17', 'Mode 18', 'Mode 19', 'Mode 20', 'Mode 21'}
%         sgtitle(join((['Heatmaps for Subject',num2str(subject),'and Emotion: ' onlineratings(trial)])))
%         filename=['export/heatmaps/people_allch/S' num2str(subject),'T' num2str(trial), extension]
%         saveas(gcf,filename)
%         close all
%         
%         MAC_plot=macmatrix(Phi_data{opt_order},Phi_cov{opt_order});
%         figure
%         macplot(MAC_plot)
%         filename=['export/heatmaps/people_allch/MAC_S' num2str(subject),'T' num2str(trial), extension]
%         saveas(gcf,filename) 
%         close all
%         
%         export_mat=[Phi_data{opt_order};zeta_data{opt_order}';fn_data{opt_order}']
%         extension='.csv'
%         filename=['export/csv_allch/S' num2str(subject),'T' num2str(trial), extension]
%         csvwrite(filename,export_mat)
%         
    end
end