set(groot, 'defaultAxesTickLabelInterpreter',"latex");
set(groot, 'defaultLegendInterpreter', "latex");
set(groot, 'defaulttextinterpreter',"latex");
colors=['b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm' 'b' 'k' 'r' 'g' 'y' 'c' 'm'];
abs_total={};
angle_total={};
fs=512;
T=1/fs;
%% Setup: DMD Heatmaps for DEAP Dataset, exports dmd heatmap and tabular version for RF

for subject = 2:2
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
    for i = 1:4
        samp_opt=[512,256,128,64];
        fact_opt=[1,2,4,8];
        n_samp=fact_opt(i);
        fs=samp_opt(i);
        Y1=s01(1:32,1:63*512);
        Y1=downsample(Y1',n_samp);
        Y1=Y1';
        
        
        extension='.png';
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
    for i = 1:4
        for ii = [46,47,48,49,50]
        abs_coef=[abs_coef;mean(mean(corrcoef(abs(Phi_out{1}(:,ii)),abs(Phi_out{i}(:,ii)))))];
        angle_coef=[angle_coef;mean(mean(corrcoef(angle(Phi_out{1}(:,ii)),angle(Phi_out{i}(:,ii)))))];
        end
    end
    abs_coef=reshape(abs_coef,[4],[]);
    angle_coef=reshape(angle_coef,[4],[]);
    abs_total{subject}=abs_coef;
    angle_total{subject}=angle_coef;
    
end



