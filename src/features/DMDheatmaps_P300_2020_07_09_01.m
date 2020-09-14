set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
filename_out=[];
export_comps=[];
export_vec=[];
fs=128; %confirmed
T=1/fs;
%% Setup: DMD Heatmaps for DEAP Dataset, exports dmd heatmap and tabular version for RF

for subject = 1:1
    if subject <= 9
        load_name1=['s0',num2str(subject),'.mat']
    else
        load_name1=['s',num2str(subject),'.mat']
    end
    subject1=subject;
    
    s01=load(load_name1);
    disp(subject);
    
    %% Loop Over Trials
    for trial = 1:1
        Y1=s01.data(trial,1:32,:);
        Y1=squeeze(Y1);
        extension='.png';
        fs=128;
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
        
%         figure('Position', [100 100 600 300]);
%         subplot(1,2,1);
%         plot(lambda, 'k.');
%         rectangle('Position', [-1 -1 2 2], 'Curvature', 1, ...
%             'EdgeColor', 'k', 'LineStyle', '--');
%         axis(1.2*[-1 1 -1 1]);
%         axis square;
%         
%         
%         subplot(1,2,2);
%         plot(omega, 'k.');
%         line([0 0], 200*[-1 1], 'Color', 'k', 'LineStyle', '--');
%         axis([-8 2 -170 +170]);
%         axis square;
%         
%         set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6 3], 'PaperPositionMode', 'manual');
%         filename2=['export/DMD/eigvals/S' num2str(subject),'T' num2str(trial), extension]
%         saveas(gcf,filename2)
%         close all
        % print('-depsc2', '-loose', ['../figures/DMD_eigenvalues.eps']);
        
        %% spectrum
        % alternate scaling of DMD modes
        Ahat = (S_r^(-1/2)) * Atilde * (S_r^(1/2));
        [What, D] = eig(Ahat);
        W_r = S_r^(1/2) * What;
        Phi = Y*V_r/S_r*W_r;
        
        f = abs(imag(omega));
        P = (diag(Phi'*Phi));
        
        % DMD spectrum
%         figure('Position', [100 100 600 300]);
%         subplot(1,2,1);
%         stem(f, P, 'k');
%         %xlim([0 150]);
%         axis square;
        
        % power spectrum
        timesteps = size(X, 2);
        srate = 1/dt;
        nelectrodes = 59;
        NFFT = 2^nextpow2(timesteps);
        f = srate/2*linspace(0, 1, NFFT/2+1);
        
%         subplot(1,2,2);
%         hold on;
%         for c = 1:nelectrodes,
%             fftp(c,:) = fft(X(c,:), NFFT);
%             plot(f, 2*abs(fftp(c,1:NFFT/2+1)), ...
%                 'Color', 0.6*[1 1 1]);
%         end;
%         plot(f, 2*abs(mean(fftp(c,1:NFFT/2+1), 1)), ...
%             'k', 'LineWidth', 2);
%         %xlim([0 150]);
%         %ylim([0 400]);
%         axis square;
%         box on;
        
%         set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6 3], 'PaperPositionMode', 'manual');
%         filename2=['export/DMD/PSD/S' num2str(subject),'T' num2str(trial), extension]
%         saveas(gcf,filename2)
        close all
        
        cutoff=length(S)/nstacks;
        Phi_phys=Phi(1:cutoff,:);
        Phi_max=max(Phi_phys);
        Phi_phys_norm=Phi_phys./Phi_max;
        Phi_phys_unique_norm=Phi_phys_norm(:,1:2:end);
        Phi_phys_unique=Phi_phys(:,1:2:end);
        t=linspace(0,timesteps/fs,timesteps);
        
        
        b = Phi\Xaug(:,1);
        for k=1:length(t)
            time_dyanmics(:,k)=(b.*exp(omega*t(k)));
        end
        Xdmd=Phi*time_dyanmics;
        
        tend=500;
        figure
        plot(t(1:tend),Y1(1,1:tend))
        hold on
        plot(t(1:tend),Y1(2,1:tend)) %select non-repeated modes (complex conj.)
        plot(t(1:tend),real(Xdmd(1,1:tend)))
        plot(t(1:tend),real(Xdmd(2,1:tend)))
%         filename2=['export/DMD/timeseries/S' num2str(subject),'T' num2str(trial), extension]
%         saveas(gcf,filename2)
%         close all
        
        figure
        set(gcf,'units','points','position',[500,-200,700,500])
        [fn,zeta]=damp(omega);
        zeta=zeta';
        fn=fn';
        fn_map=fn/max(fn);
        fn_map=fn_map(:,1:2:end);
        zeta_map=zeta/max(zeta);
        zeta_map=zeta_map(:,1:2:end);
        h_indmap=heatmap([real(Phi_phys_unique_norm);imag(Phi_phys_unique_norm);zeta_map;fn_map],'CellLabelColor','none');
        h_indmap.Colormap=parula;
        h_indmap.ColorbarVisible=0;        
        %sgtitle(join((['Heatmaps (Averages of A) for Emotion: ' onlineratings(trial)])))
        filename2=['export/DMD_tab/robots/S' num2str(subject),'T' num2str(trial), extension]
        saveas(gcf,filename2)
        close all
        
        subject_vec=repelem(subject,size([Phi_phys(:,1:2:end);zeta(:,1:2:end);fn(:,1:2:end)],1));
        trial_vec=repelem(trial,size([Phi_phys(:,1:2:end);zeta(:,1:2:end);fn(:,1:2:end)],1));
        
        Phi_phys_exp=Phi_phys(:,1:2:end);
        Phi_phys_exp=Phi_phys_exp(:);
        Phi_phys_exp_re=real(Phi_phys_exp);
        Phi_phys_exp_im=imag(Phi_phys_exp);
        subject_exp=repelem(subject,(size(Phi_phys,1))*r/2,1);
        trial_exp=repelem(trial,(size(Phi_phys,1))*r/2,1);
        chann_count=(1:size(Phi_phys,1))';
        chann_exp=repmat(chann_count,r/2,1);
        
        mode_count=(1:size(Phi_phys,2)/2)';
        mode_exp=repelem(mode_count,size(Phi_phys,1),1);
        
        fn_count=fn(:,1:2:end);
        fn_exp=repelem(fn_count',size(Phi_phys,1),1);
        
        zeta_count=zeta(:,1:2:end);
        zeta_exp=repelem(zeta_count',size(Phi_phys,1),1);
        
        export_vec=[export_vec;[subject_exp(1:r/2),trial_exp(1:r/2),mode_count,real(Phi_phys_unique)',imag(Phi_phys_unique)',fn_count',zeta_count']];
        export_comps=[export_comps;[subject_exp,trial_exp,Phi_phys_exp_re,Phi_phys_exp_im,chann_exp,mode_exp,fn_exp,zeta_exp]];

    end
end
writematrix(export_vec,'dmd_deap_100modes_vecs.csv')
writematrix(export_comps,'dmd_deap_100modes_comps.csv')