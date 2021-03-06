set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
filename_out=[];
fs=160;
T=1/fs;
%% Setup: DMD Heatmaps for DEAP Dataset
for subject = 1:109
    load_name_mat=[];
    if subject <= 9
        for trial = 1:9            
            load_name = ['S00',num2str(subject),'R0',num2str(trial),'.h5'];
            load_name_mat=[load_name_mat;load_name];
        end
        for trial = 10:14
            load_name = ['S00',num2str(subject),'R',num2str(trial),'.h5'];
            load_name_mat=[load_name_mat;load_name];
        end
    elseif subject <=99 && subject >=10
        for trial = 1:9            
            load_name = ['S0',num2str(subject),'R0',num2str(trial),'.h5'];
            load_name_mat=[load_name_mat;load_name];
        end
        for trial = 10:14
            load_name = ['S0',num2str(subject),'R',num2str(trial),'.h5'];
            load_name_mat=[load_name_mat;load_name];
        end
    elseif subject >= 100
        for trial = 1:9            
            load_name = ['S',num2str(subject),'R0',num2str(trial),'.h5']
            load_name_mat=[load_name_mat;load_name];
        end
        for trial = 10:14
            load_name = ['S',num2str(subject),'R',num2str(trial),'.h5']
            load_name_mat=[load_name_mat;load_name];
        end
    end
    
    %% Loop Over Trials
    for trial = 1:14
        s01 = h5read(load_name_mat(trial,:),['/df/block0_values']); %
        fix=length(s01)/4;
        Y_split{1}=s01(1:32,1:fix*1);
        Y_split{2} = s01(1:32,fix:fix*2);
        Y_split{3} = s01(1:32,fix*2:fix*3);
        extension='.png';
        fs=160;
        dt=1/fs;
        
        parfor i = 1:3
            Y_curr=Y_split{i};
            r = 80; % number of modes, remember DMD generates complex conjugates
            nstacks = 10; % number of stacks
        
            % construct the augmented, shift-stacked data matrices
            Xaug = [];
            for st = 1:nstacks,
                Xaug = [Xaug; Y_curr(:, st:end-nstacks+st)]; %
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
            filename2=['export/DMD_phys/eigvals/S' num2str(subject),'T' num2str(trial*i), extension]
            saveas(gcf,filename2)
            close all
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
            figure('Position', [100 100 600 300]);
            subplot(1,2,1);
            stem(f, P, 'k');
            %xlim([0 150]);
            axis square;

            % power spectrum
            timesteps = size(X, 2);
            srate = 1/dt;
            nelectrodes = 59;
            NFFT = 2^nextpow2(timesteps);
            f = srate/2*linspace(0, 1, NFFT/2+1);

            subplot(1,2,2);
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

            set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6 3], 'PaperPositionMode', 'manual');
            filename2=['export/DMD_phys/PSD/S' num2str(subject),'T' num2str(trial*i), extension]
            saveas(gcf,filename2)
            close all

            cutoff=length(S)/nstacks;
            Phi_phys=Phi(1:cutoff,:);
            Phi_max=max(Phi_phys);
            Phi_phys_norm=Phi_phys./Phi_max;
            Phi_phys_unique=Phi_phys_norm(:,1:2:end);
            t=linspace(0,timesteps/fs,timesteps);


            b = Phi\Xaug(:,1);
    %         for k=1:length(t)
    %             time_dyanmics(:,k)=(b.*exp(omega*t(k)));
    %         end
    %         Xdmd=Phi*time_dyanmics;

    %         tend=500;
    %         figure
    %         plot(t(1:tend),Y1(1,1:tend))
    %         hold on
    %         plot(t(1:tend),Y1(2,1:tend)) %select non-repeated modes (complex conj.)
    %         plot(t(1:tend),real(Xdmd(1,1:tend)))
    %         plot(t(1:tend),real(Xdmd(2,1:tend)))
    %         filename2=['export/DMD/timeseries/S' num2str(subject),'T' num2str(trial), extension]
    %         saveas(gcf,filename2)
    %         close all

            figure
            set(gcf,'units','points','position',[500,-200,700,500])
            fn_map=abs(transpose(omega/max(omega)));
            fn_map=fn_map(:,1:2:end);
            h_indmap=heatmap([real(Phi_phys_unique);imag(Phi_phys_unique);fn_map],'CellLabelColor','none');
            h_indmap.Colormap=parula;
            h_indmap.ColorbarVisible=0;        
            %sgtitle(join((['Heatmaps (Averages of A) for Emotion: ' onlineratings(trial)])))
            filename2=['export/DMD_phys/robots/S' num2str(subject),'T' num2str(trial*i), extension]
            saveas(gcf,filename2)
            close all

            %heatmaps for people
    %         figure
    %         set(gcf,'units','points','position',[500,-200,700,500])
    %         
    %         h_indmap=heatmap([real(Phi_phys_unique);imag(Phi_phys_unique);fn_map],'CellLabelColor','none');
    %         h_indmap.Colormap=parula;
    %         h_indmap.ColorbarVisible=0;
    %         %h_indmap.XDisplayLabels={'Mode 1', 'Mode 2', 'Mode 3', 'Mode 4', 'Mode 5', 'Mode 6', 'Mode 7', 'Mode 8', 'Mode 9', 'Mode 10', 'Mode 11', 'Mode 12', 'Mode 13', 'Mode 14', 'Mode 15', 'Mode 16', 'Mode 17', 'Mode 18', 'Mode 19', 'Mode 20'}
    %         sgtitle(join((['Heatmaps for Subject',num2str(subject),'and Emotion: ' onlineratings(trial)])))
    %         filename=['export/DMD/people/S' num2str(subject),'T' num2str(trial), extension]
    %         saveas(gcf,filename)
    %         close all

        end

    end
    
    
end