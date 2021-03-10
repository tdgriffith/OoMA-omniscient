set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
filename_out=[];
fs=160;
T=1/fs;
extension='.png';
%% Setup: DMD Heatmaps for Physio Dataset
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
    count=0;
    for trial = 3:14
        s01 = h5read(load_name_mat(trial,:),['/df/block0_values']); %
        s01=reshape(s01,64,[],4);
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
        
        for i = 1:size(s01,3)
            Y_curr=s01(:,:,i);
            count=count+1;
            % parameters:
            r = 200; % number of modes, remember DMD generates complex conjugates
            nstacks = 10; % number of stacks
            
            % construct the augmented, shift-stacked data matrices
            Xaug = [];
            for st = 1:nstacks,
                Xaug = [Xaug; Y_curr(:, st:end-nstacks+st)];%#ok<AGROW>
            end;
            X = Xaug(:, 1:end-1);
            Y = Xaug(:, 2:end);
            
            % SVD and truncate to first r modes
            [U, S, V] = svd(X, 'econ');
            U_r = U(:, 1:r);
            S_r = S(1:r, 1:r);
            V_r = V(:, 1:r);
            
            %Singular Value Plots
            %%

            
            % DMD modes
            Atilde = U_r'*Y*V_r/S_r;
            [W_r, D] = eig(Atilde);
            Phi = Y*V_r/S_r*W_r;
            
            % DMD eigenvalues
            lambda = diag(D);
            omega = log(lambda)/dt/2/pi;
            Omega=diag(omega);
            
            %% eigenvalue

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
            Phi_phys_unique=Phi_phys_norm(:,1:2:end);
            t=linspace(0,timesteps/fs,timesteps);
            
            

            
            tend=500;
            
            figure
            set(gcf,'units','points','position',[500,-200,700,500])
            fn_map=abs(transpose(omega/max(omega)));
            fn_map=fn_map(:,1:2:end);
            h_indmap=heatmap([real(Phi_phys_unique);imag(Phi_phys_unique);fn_map],'CellLabelColor','none');
            h_indmap.Colormap=parula;
            h_indmap.ColorbarVisible=0;
            filename2=['export/DMD_phys_3/S' num2str(subject),'_trial',num2str(trial),'Task' num2str(task),'_i',num2str(i), extension]
            saveas(gcf,filename2)
            close all


        end

    end
    
    
    
end