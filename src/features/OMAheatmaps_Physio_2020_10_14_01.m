set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
filename_out=[];
export_mat=[];
fs=160;
T=1/fs;
%% Setup: DMD Heatmaps for Physio Dataset
for subject = 1:40 %109
    load_name_mat=[];
    for trial = 1:9
        load_name = ['S',num2str(subject),'R0',num2str(trial),'.h5']
        load_name_mat=[load_name_mat;load_name];
    end
    for trial = 10:14
        load_name = ['S',num2str(subject),'R',num2str(trial),'.h5']
        load_name_mat=[load_name_mat;load_name];
    end

    
    %% Loop Over Trials
    count=0;
    for trial = 3:14 %14
        s01 = h5read(load_name_mat(trial,:),['/df/block0_values']); %
        data=zeros(64,3200,7);
        for n=1:size(s01,1)
            data(n,:,:)=buffer(s01(n,:),20*160,2*160);
        end
        if trial == 3||7||11
            task=1
        elseif trial == 4||6||12
            task=2
        elseif trial == 5||9||13
            task=3
        elseif trial == 6||10||14
            task=4
        end
        extension='.png';
        fs=160;
        dt=1/fs;
        
        for i = 1:size(data,3)
            Y_curr=data(1:64,:,i);
            Y_filt=bandpass(Y_curr',[4 45],fs);
            Y_filt=Y_filt';
            %count=count+1;
            order = 32;
            s = 2*order;
            opt_order=32;
            [A_data,C_data,G_data,R0_data] = ssicov(Y_filt,order,s);
            err = [0.01,0.05,0.98];

            
            [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,dt);
            Phi_data=Phi_data';
            norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
            Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest
            %[IDs_cov] = plotstab(A_data,C_data,Y_filt,dt,[],err);
            
            
            
            %heatmaps for robuts
            figure
            set(gcf,'units','points','position',[500,-200,700,500])
            zeta_map=transpose(zeta_data{opt_order}/max(zeta_data{opt_order}));
            fn_map=transpose(fn_data{opt_order}/max(fn_data{opt_order}));
            h_indmap=heatmap([real(Phi_data{opt_order});imag(Phi_data{opt_order});zeta_map;fn_map]);
            h_indmap.Colormap=parula;
            h_indmap.ColorbarVisible=0;
            
            %sgtitle(join((['Heatmaps (Averages of A) for Emotion: ' onlineratings(trial)])))
            filename2=['export/eegmmi/S' num2str(subject),'T' num2str(trial),'W',num2str(i), extension]
            saveas(gcf,filename2)
            close all
            export_vec=[subject,trial,reshape([real(Phi_data{opt_order}(:,1:15));imag(Phi_data{opt_order}(:,1:15));zeta_data{opt_order}(1:15).';fn_data{opt_order}(1:15).'],1,[])];
            export_mat=[export_mat;export_vec];


        end

    end
    
    
    
end
writematrix(export_mat,'OMA_EEGmmi_window_nohead.csv')
%writematrix(export_comps,'dmd_deap_100modes_comps.csv')
header={};
for mode = 1:size(Phi_data{opt_order}(:,1:15),2)
    for comp = 1:size(Phi_data{opt_order}(:,1:15),1)
        header = [header,{['Real_Comp',num2str(comp),'_Mode',num2str(mode)]}];
    end
end

for mode = 1:size(Phi_data{opt_order}(:,1:15),2)
    for comp = 1:size(Phi_data{opt_order}(:,1:15),1)
        header = [header,{['Imag_Comp',num2str(comp),'_Mode',num2str(mode)]}];
    end
end

for mode = 1:size(Phi_data{opt_order}(:,1:15),2)
    header = [header,{['zeta',num2str(mode)]}];
end

for mode = 1:size(Phi_data{opt_order}(:,1:15),2)
    header = [header,{['fn',num2str(mode)]}];
end

header = [{'Subject'},{'Trial'},header];

out = [header;num2cell(export_mat)];
writecell(out,'OMA_EEGmmi_window_head.csv')