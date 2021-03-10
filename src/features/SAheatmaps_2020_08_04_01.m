set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
filename_out=[];
export_comps=[];
export_vec=[];
fs=128;
T=1/fs;
%% Setup: DMD Heatmaps for DEAP Dataset, exports dmd heatmap and tabular version for RF
channel_cov=readtable('deap_imag_conv.csv');
for subject = 1:1
    if subject <= 9
        load_name1=['s0',num2str(subject),'.mat'];
    else
        load_name1=['s',num2str(subject),'.mat'];
    end
    subject1=subject;
    
    s01=load(load_name1);
    disp(subject);
    
    %% Loop Over Trials
    for trial = 1:1
        %load data
        Y1=s01.data(trial,1:32,:);
        Y1=squeeze(Y1);
        Y1 = Y1';
        extension='.png';
        fs=128;
        dt=1/fs;
        %t=[dt:dt:(size(Y1,1)*dt)];
        
        %calculate derivative
        dy = gradient(Y1,dt);
        dy2=reshape(dy.',1,[])';
        
        %Define library
        n=size(Y1,2);
        basis=gen_SA_basis_noj(n);
        Theta_cust=[];
        for i = 1:length(basis)
            Y_sigma = Y1*(basis{i});
            Y_sigma = reshape(Y_sigma.',1,[])';
            Theta_cust = [Theta_cust,Y_sigma];
        end
        
        Xi = [real(Theta_cust);imag(Theta_cust)]\[real(dy2);imag(dy2)]; %Really not sure here. Something messed up with imag why backwards???
        
        basis=gen_SA_basis_noj2(n);
        A_solve=zeros(n,n);
        for i = 1:length(basis)
            A_solve=A_solve+(Xi(i)*basis{i});
        end   
        
        H_solve=A_solve/(-1j);
        
        [Phi,d]=eig(H_solve);
        Phi_max=max(Phi);
        Phi_norm=Phi./Phi_max;
        
        figure
        set(gcf,'units','points','position',[500,-200,700,500])
        [fn,zeta]=damp(A_solve);
        zeta=zeta';
        fn=fn';
        fn_map=fn/max(fn);
        %zeta_map=zeta/max(zeta);
        %h_indmap=heatmap([real(Phi_norm);imag(Phi_norm);zeta_map;fn_map],'CellLabelColor','none');
        h_indmap=heatmap([real(Phi_norm);imag(Phi_norm);fn_map],'CellLabelColor','none');
        h_indmap.Colormap=parula;
        h_indmap.ColorbarVisible=0;        
        %sgtitle(join((['Heatmaps (Averages of A) for Emotion: ' onlineratings(trial)])))
        filename2=['export/SA_mats/robots/S' num2str(subject),'T' num2str(trial), extension]
        saveas(gcf,filename2)
        close all
        
        %HERE HERE
        subject_vec=repelem(subject,size([Phi;fn],1));
        trial_vec=repelem(trial,size([Phi;fn],1));
        
        Phi_phys_exp=Phi;
        Phi_phys_exp=Phi_phys_exp(:);
        Phi_phys_exp_re=real(Phi_phys_exp);
        Phi_phys_exp_im=imag(Phi_phys_exp);
        subject_exp=repelem(subject,(size(Phi,1)),1);
        trial_exp=repelem(trial,(size(Phi,1)),1);
        chann_count=(1:size(Phi,1))';
        chann_exp=repmat(chann_count,1);
        
        mode_count=(1:size(Phi,2))';
        mode_exp=repelem(mode_count,size(Phi,1),1);
        
        fn_count=fn();
        fn_exp=repelem(fn_count',size(Phi,1),1);
        
   
        
      
        
        
        %export_vec=[export_vec;reshape([subject_exp(1:r/2),trial_exp(1:r/2),mode_count,real(Phi_phys_unique)',imag(Phi_phys_unique)',fn_count',zeta_count'],1,[])];
        export_vec=[export_vec;[subject,trial,reshape([real(Phi)',imag(Phi)',fn_count'],1,[])]];
        %export_comps=[export_comps;[subject_exp,trial_exp,Phi_phys_exp_re,Phi_phys_exp_im,chann_exp,mode_exp,fn_exp,zeta_exp]];

    end
end
writematrix(export_vec,'2020-08-20_SA_deap_modes_nohead.csv')

header={};
for mode = 1:size(Phi,2)
    for comp = 1:size(Phi,1)
        header = [header,{['Real_Comp',num2str(comp),'_Mode',num2str(mode)]}];
    end
end

for mode = 1:size(Phi,2)
    for comp = 1:size(Phi,1)
        header = [header,{['Imag_Comp',num2str(comp),'_Mode',num2str(mode)]}];
    end
end

for mode = 1:size(Phi,2)
    header = [header,{['fn',num2str(mode)]}];
end

header = [{'Subject'},{'Trial'},header];

out = [header;num2cell(export_vec)];
writecell(out,'2020-08-20_SA_deap_modes_trials.csv')


