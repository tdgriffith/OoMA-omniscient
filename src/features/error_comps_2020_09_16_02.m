set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
abs_total={};
angle_total={};
fs=512;
T=1/fs;
%% Setup: DMD Heatmaps for DEAP Dataset, exports dmd heatmap and tabular version for RF

parfor subject = 1:32
    if subject <= 9
        load_name1=['s0',num2str(subject),'_512.h5']
    else
        load_name1=['s',num2str(subject),'_512.h5']
    end
    
    s01=h5read(load_name1,'/df_512/block0_values');
    disp(subject);
    
    %% Loop Over Trials
    Phi_out={};
    fn_out={};
    for i = 1:8
        snr_opt=[800,100,40,20,10,5,2,1];
        noise_int=snr_opt(i)
        
        Y1=s01(1:32,1:63*512);
        Y1=awgn(Y1',noise_int,'measured');
        Y1=Y1';
        if noise_int >100
            Y1=s01(1:32,1:63*512);
            display('Skipped')
        end
        extension='.png';
        fs=512;
        dt=1/fs;
        
        order = 35;
        s = 2*order;
        opt_order=35;
        % OMA-data
        [A_data,C_data,G_data,R0_data] = ssidata(Y1,order,s);
        err = [0.01,0.05,0.98];
        
        
        [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T);
        Phi_data=Phi_data';
       
        Phi_out{i}=Phi_data{opt_order};
        
        
        fn_out{i}=fn_data{opt_order};

    end
    abs_coef=[];
    angle_coef=[];
    for i = 1:8
        for ii = [1,5,10,15]
        abs_coef=[abs_coef;mean(mean(corrcoef(abs(Phi_out{1}(:,ii)),abs(Phi_out{i}(:,ii)))))];
        angle_coef=[angle_coef;mean(mean(corrcoef(angle(Phi_out{1}(:,ii)),angle(Phi_out{i}(:,ii)))))];
        end
    end
    abs_coef=reshape(abs_coef,[4],[]);
    angle_coef=reshape(angle_coef,[4],[]);
    abs_total{subject}=abs_coef;
    angle_total{subject}=angle_coef;
    
    
end
abs_res = cat(3,abs_total{:});
abs_res_mean = mean(abs_res,3);

angle_res = cat(3,angle_total{:});
angle_res_mean = mean(angle_res,3);


x=[0,1/80,1/40,1/20,1/10,1/5,1/2,1];
figure
subplot(1,2,1)
semilogx(x,abs_res_mean)
ylim([0,1])
%set(gca, 'XDir','reverse')
subplot(1,2,2)
semilogx(x,angle_res_mean)
ylim([0,1])
%set(gca, 'XDir','reverse')
% writematrix(export_vec,'dmd_deap_100modes_trials3_nohead.csv')
% %writematrix(export_comps,'dmd_deap_100modes_comps.csv')
% header={};
% for mode = 1:50
%     for comp = 1:11
%         header = [header,{['Real_Comp',num2str(comp),'_Mode',num2str(mode)]}];
%     end
% end
% 
% for mode = 1:50
%     for comp = 1:11
%         header = [header,{['Imag_Comp',num2str(comp),'_Mode',num2str(mode)]}];
%     end
% end
% 
% for mode = 1:50
%     header = [header,{['fn',num2str(mode)]}];
% end
% 
% for mode = 1:50
%     header = [header,{['zeta',num2str(mode)]}];
% end
% 
% header = [{'Subject'},{'Trial'},header];
% 
% out = [header;num2cell(export_vec)];
% writecell(out,'dmd_deap_100modes_trials3.csv')


