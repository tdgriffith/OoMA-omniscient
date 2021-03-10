set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
%onlineratings=gen_onlineratings();
filename_out=[];
export_vec=[];
%% Setup: Heatmaps for Robots NORMALIZED OMA DATA
for subject = 1:32 %parfor
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
    s02=load(load_name2)
    
    %mkdir(basename)

%% Loop Over Trials
    for trial = 1:40
        % Trial Data
        %Y1=s01.data(trial,[3,4,31,27,28],:);
        Y1=s01.data(trial,1:15,:);
        Y1=squeeze(Y1);
        %Y2=s01.data(trial+1,[3,4,31,27,28],:);
        %Y2=squeeze(Y2);
        t = linspace(0,63,length(Y1));
        fs=128;
        T=1/fs;
        dt=T;
        extension='.png';

        % OMA Covariance Algorithm
        order = 40;
        s = 2*order;
        opt_order=30;
        %[A_cov,C_cov,G_cov,R0_cov,S_cov] = ssicov(Y1,order,s);

        %eig(A_cov{opt_order})
        err = [0.01,0.05,0.98];
        

        [A_data,C_data,G_data,R0_data] = ssidata(Y1,order,s);
        err = [0.01,0.05,0.98];

        
        [fn_data,zeta_data,Phi_data] = modalparams(A_data,C_data,T);
        Phi_data=Phi_data';
        norm_data=abs(max(Phi_data{opt_order}, [], 'all')); %determine max value
        Phi_data{opt_order}=Phi_data{opt_order}/norm_data; %normalize the one of interest
        
        
        export_vec=[export_vec;[subject,trial,reshape([real(eig(A_data{opt_order}))',imag(eig(A_data{opt_order}))'],1,[])]];
        lambda=eig(A_data{opt_order});
        omega=log(lambda)/dt/2/pi;
        figure('Position', [100 100 600 300]);
        subplot(1,2,1);
        plot(lambda, 'k.');
        rectangle('Position', [-1 -1 2 2], 'Curvature', 1, ...
            'EdgeColor', 'k', 'LineStyle', '--');
        axis(1.2*[-1 1 -1 1]);
        axis square;
        
        
        subplot(1,2,2);
        plot(omega, 'k.');
        line([0 0], 200*[-1 1], 'Color', 'k', 'LineStyle', '--');
        axis([-8 2 -170 +170]);
        axis square;
        set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6 3], 'PaperPositionMode', 'manual');
        mkdir('export/OMA_lambda')
        extension='.png'
        filename2=['/mnt/tris_files/export/OMA_lambda/S' num2str(subject),'T' num2str(trial), extension]
        saveas(gcf,filename2)
        close all        
    end
end
writematrix(export_vec,'/mnt/tris_files/export/OMA_deap_30order_eigvals_nohead.csv')
%writematrix(export_comps,'dmd_deap_100modes_comps.csv')
header={};
for mode = 1:opt_order
        header = [header,{['Real_Comp',num2str(mode),'_Mode',num2str(mode)]}];
end

for mode = 1:opt_order
        header = [header,{['Imag_Comp',num2str(mode),'_Mode',num2str(mode)]}];
end



header = [{'Subject'},{'Trial'},header];
out = [header;num2cell(export_vec)];
writecell(out,'/mnt/tris_files/export/OMA_deap_30order_eigvals_head.csv')