set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
onlineratings=gen_onlineratings();
filename_out=[];
fs=512;
T=1/fs;
%% Setup: Heatmaps for Robots LH RH Data set
% Hohyun Cho, Minkyu Ahn, Sangtae Ahn, Moonyoung Kwon, Sung Chan Jun, 
% EEG datasets for motor imagery braincomputer interface, GigaScience, Volume 6, Issue 7, 
% July 2017, gix034, https://doi.org/10.1093/gigascience/gix034

parfor subject = 1:52
    if subject <= 9
        load_name1=['s0',num2str(subject),'.mat']
    else
        load_name1=['s',num2str(subject),'.mat']
    end
    subject1=subject;

    s01=load(load_name1);
    s01=s01.eeg.rest;
    s01=s01(1:64,1:30720);
    basename=['S' num2str(subject1)];
    extension='.png';
    s01=detrend(s01');
    s01=bandstop(s01,[59 61],fs);
    s01=bandpass(s01,[3, 50],fs);
    
    s01=s01';
    s01=reshape(s01,64,7680,4);

%% Loop Over Trials
    for trial = 1:4
        Y1=s01(:,:,trial);
        Y1=squeeze(Y1);
        extension='.png';

        % OMA Covariance Algorithm
        order = 25;
        s = 2*order;
        opt_order=20;
        [A_cov,C_cov,G_cov,R0_cov,S_cov] = ssicov(Y1,order,s);
        err = [0.01,0.05,0.98];
        %[IDs_cov] = plotstab(A_cov,C_cov,Y1,T,[],err);
        [fn_cov,zeta_cov,Phi_cov] = modalparams(A_cov,C_cov,T);
        norm_cov=abs(max(Phi_cov{opt_order}, [], 'all'));
        Phi_cov{opt_order}=Phi_cov{opt_order}/norm_cov;
      


        % OMA-data
        [A_data,C_data,G_data,R0_data] = ssidata(Y1,order,s);
        err = [0.01,0.05,0.98];
        %[IDs_cov] = plotstab(A_data,C_data,Y1,T,[],err);

        [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T);
        Phi_data=Phi_data';
        norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
        Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest



        
        %heatmaps for robuts
        figure
        set(gcf,'units','points','position',[500,-200,700,500])
        zeta_map=transpose(zeta_data{opt_order}/max(zeta_data{opt_order}));
        fn_map=transpose(fn_data{opt_order}/max(fn_data{opt_order}));
        h_indmap=heatmap([real(Phi_data{opt_order});imag(Phi_data{opt_order});zeta_map;fn_map]);
        h_indmap.Colormap=parula;
        h_indmap.ColorbarVisible=0;
        
        %sgtitle(join((['Heatmaps (Averages of A) for Emotion: ' onlineratings(trial)])))
        filename2=['export/heatmaps/resting2/S' num2str(subject),'T' num2str(trial), extension]
        saveas(gcf,filename2)
        close all
        
        %heatmaps for people
        figure
        set(gcf,'units','points','position',[500,-200,700,500])
        
        h_indmap=heatmap([real(Phi_data{opt_order});imag(Phi_data{opt_order});zeta_map;fn_map]);
        h_indmap.Colormap=parula;
        h_indmap.ColorbarVisible=0;
        %h_indmap.XDisplayLabels={'Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Mode 5', 'Mode 6', 'Mode 7', 'Mode 8', 'Mode 9', 'Mode 10', 'Mode 11', 'Mode 12', 'Mode 13', 'Mode 14', 'Mode 15', 'Mode 16', 'Mode 17', 'Mode 18', 'Mode 19', 'Mode 20'}
        sgtitle(join((['Heatmaps for Subject',num2str(subject),'and Emotion: ' onlineratings(trial)])))
        filename=['export/heatmaps/resting2/S' num2str(subject),'T' num2str(trial), extension]
        saveas(gcf,filename)
        close all
        
        MAC_plot=macmatrix(Phi_data{opt_order},Phi_cov{opt_order});
        figure
        macplot(MAC_plot)
        filename=['export/heatmaps/resting2/MAC_S' num2str(subject),'T' num2str(trial), extension]
        saveas(gcf,filename) 
        close all
        
        output_mat=[real(Phi_data{opt_order});imag(Phi_data{opt_order});zeta_data{opt_order}';fn_data{opt_order}'];
        output_vec=output_mat(:);
        extension='.csv'
        filename=['export/csv_resting/S' num2str(subject),'T' num2str(trial), extension]
        csvwrite(filename,output_vec)
        
    end
end