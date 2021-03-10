%% TODO OVERLAPPING WINDOWS, EXPORT ALL AND LEAVE NANS if fewer modes extracted; fix tabular header

set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm''b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
onlineratings=gen_onlineratings();
filename_out=[];
export_mat=[];
fs=128;
T=1/fs;
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
    data=s01.data(:,1:32,:);
    data2=data;
%     data2=reshape(data,160,32,[]);


    
    %mkdir(basename)

%% Loop Over Trials
    for trial = 1:size(data2,1)
        % Trial Data
        %Y1=s01.data(trial,[3,4,31,27,28],:);
        Y1=data2(trial,:,:);
        Y1=squeeze(Y1);
        %t = linspace(0,size(data,1)/fs,length(Y1));
        extension='.png';

        % OMA Covariance Algorithm
        order = 15;
        s = 2*order;
        opt_order=15;
       

        % OMA-data
        [A_data,C_data,G_data,R0_data] = ssidata(Y1,order,s);
        err = [0.01,0.05,0.98];
        %[IDs_cov] = plotstab(A_data,C_data,Y1,T,[],err);

        
        [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T);
        Phi_data=Phi_data';
        norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
        Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest


        


        
        %heatmaps for robuts
%         figure
%         set(gcf,'units','points','position',[500,-200,700,500])
%         zeta_map=transpose(zeta_data{opt_order}/max(zeta_data{opt_order}));
%         fn_map=transpose(fn_data{opt_order}/max(fn_data{opt_order}));
%         h_indmap=heatmap([real(Phi_data{opt_order});imag(Phi_data{opt_order});zeta_map;fn_map]);
%         h_indmap.Colormap=parula;
%         h_indmap.ColorbarVisible=0;
%         
%         %sgtitle(join((['Heatmaps (Averages of A) for Emotion: ' onlineratings(trial)])))
%         filename2=['export/OMA_window/robots/S' num2str(subject),'T' num2str(trial), extension]
%         saveas(gcf,filename2)
%         close all
        

%         
        export_vec=[repelem(trial,size(Phi_data{opt_order},2));repelem(trial,size(Phi_data{opt_order},2));linspace(1,size(Phi_data{opt_order},2),size(Phi_data{opt_order},2));real(Phi_data{opt_order});imag(Phi_data{opt_order});zeta_data{opt_order}.';fn_data{opt_order}.'].';
        export_mat=[export_mat;export_vec];
%         extension='.csv'
%         filename=['export/csv_allch/S' num2str(subject),'T' num2str(trial), extension]
%         csvwrite(filename,export_mat)
%         
    end
end
writematrix(export_mat,'OMA_deap_window_nohead2.csv')
%writematrix(export_comps,'dmd_deap_100modes_comps.csv')
header={};

for comp = 1:size(Phi_data{opt_order},1)
    header = [header,{['Real_Comp',num2str(comp)]}];
end



for comp = 1:size(Phi_data{opt_order},1)
    header = [header,{['Imag_Comp',num2str(comp)]}];
end


% for mode = 1:size(Phi_data{opt_order}(:,1:7),2)
%     header = [header,{['zeta',num2str(mode)]}];
% end
% 
% for mode = 1:size(Phi_data{opt_order}(:,1:7),2)
%     header = [header,{['fn',num2str(mode)]}];
% end

header = [{'Subject'},{'Trial'},{'Mode No.'},header,{'zeta'},{'fn'}];

out = [header;num2cell(export_mat)];
writecell(out,'OMA_deap_window_head2.csv')