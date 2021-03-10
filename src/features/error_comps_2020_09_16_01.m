set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
abs_total={};
angle_total={};
fs=512;
T=1/fs;
%% Setup: DMD Heatmaps for DEAP Dataset, exports dmd heatmap and tabular version for RF

for subject = 1:32
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
        
        r = 100; % number of modes, remember DMD generates complex conjugates
        nstacks = 10; % number of stacks
        
        % construct the augmented, shift-stacked data matrices
        Xaug = [];
        for st = 1:nstacks,
            Xaug = [Xaug; Y1(:, st:end-nstacks+st)]; %
        end;
        X = Xaug(:, 1:end-1);
        Y = Xaug(:, 2:end);
        
        % SVD and truncate to first r modes
        [U, S, V] = svd(X, 'econ');
        U_r = U(:, 1:r);
        S_r = S(1:r, 1:r);
        V_r = V(:, 1:r);
        
        % DMD modes
        Atilde = U_r'*Y*V_r/S_r;
        [W_r, D] = eig(Atilde);
        Phi = Y*V_r/S_r*W_r;
        
        % DMD eigenvalues
        lambda = diag(D);
        omega = log(lambda)/dt/2/pi;
        Omega=diag(omega);
        
        
        %% spectrum
        % alternate scaling of DMD modes
        Ahat = (S_r^(-1/2)) * Atilde * (S_r^(1/2));
        [What, D] = eig(Ahat);
        W_r = S_r^(1/2) * What;
        Phi = Y*V_r/S_r*W_r;
        
        f = abs(imag(omega));
        P = (diag(Phi'*Phi));
        

        
        % power spectrum
        timesteps = size(X, 2);
        srate = 1/dt;
        nelectrodes = 59;
        NFFT = 2^nextpow2(timesteps);
        f = srate/2*linspace(0, 1, NFFT/2+1);
        

        
        cutoff=length(S)/nstacks;
        Phi_phys=Phi(1:cutoff,:);
        Phi_max=max(Phi_phys);
        Phi_phys_norm=Phi_phys./Phi_max;
        Phi_phys_unique_norm=Phi_phys_norm(:,1:2:end);
        Phi_phys_unique=Phi_phys(:,1:2:end);
        Phi_out{i}=Phi_phys_unique;
        t=linspace(0,timesteps/fs,timesteps);
        
        
        b = Phi\Xaug(:,1);
        for k=1:length(t)
            time_dyanmics(:,k)=(b.*exp(omega*t(k)));
        end
        Xdmd=Phi*time_dyanmics;
        
        tend=500;
        
        [fn,zeta]=damp(omega);
        zeta=zeta';
        fn=fn';
        fn=fn(:,1:2:end);
        fn_out{i}=fn;

    end
    abs_coef=[];
    angle_coef=[];
    for i = 1:8
        for ii = [20,30,35,40]
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


x=[0,1/100,1/40,1/20,1/10,1/5,1/2,1];
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


