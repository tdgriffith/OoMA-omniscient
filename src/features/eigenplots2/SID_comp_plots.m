%% Single Subject, All Trials, Compare SID Methods
%% Setup 
loyolagreen=1/255*[0,104,87];
for subject = 1:32
    if subject <= 9
        load_name=['s0',num2str(subject),'.mat']
    else
        load_name=['s',num2str(subject),'.mat']
    end        
    s02=load(load_name);

    basename='Algorithm Comparison for Subject';
    extension='.png';

    figure
    set(gcf,'units','points','position',[100,100,700,700])
    filename=[basename, num2str(subject), extension]
    %% Loop over Trials
    for trial = 1:40

        Y=s02.data(trial,1:6,:);
        Y=squeeze(Y);
        t = linspace(0,63,length(Y));
        fs=128;
        T=1/fs;


        % TS and PS Plots, don't need rn.
        %[Syy,freqs] = pwelch(Y',[],[],[],fs); % obtain estimates of the output power spectrums
        % clf
        % subplot(2,1,1)
        % plot(t,Y)
        % xlabel('Time (s)')
        % ylabel('EEG Output (uV)')
        % axis tight
        % subplot(2,1,2)
        % plot(freqs,10*log10(Syy))
        % xlabel('Frequency (Hz)')
        % ylabel('PSD (dB)')
        % axis tight

        % OMA Algorithm (Covar based)
        order = 25;
        s = 2*order;
        opt_order=12; %purely just cause here, seems to have the most modes
        [A_cov,C_cov,G_cov,R0_cov] = ssicov(Y,order,s);
        eig(A_cov{opt_order})
        err = [0.01,0.05,0.98];
        %[IDs_cov] = plotstab(A_cov,C_cov,Y,T,[],err);
        %[fn_cov,zeta_cov,Phi_cov] = modalparams(A_cov,C_cov,T);

        real_part_OMA=real(eig(A_cov{opt_order}));
        imag_part_OMA=imag(eig(A_cov{opt_order}));

        %n4sid
        data = iddata(Y',[],T);
        sys = n4sid(data,12);
        [A_n4,~,C_n4,~] = ssdata(sys);
        real_part_n4=real(eig(A_n4));
        imag_part_n4=imag(eig(A_n4));    

        % OMA-data
        [A_data,C_data,G_data,R0_data] = ssidata(Y,order,s);
        real_part_data=real(eig(A_data{opt_order}));
        imag_part_data=imag(eig(A_data{opt_order}));
        
        %NeXT-ERA
        [NeXT]=NExTFERA(Y,6,2000,4,0.1,fs,800,200,12,10,1);
        real_part_NeXT=real(eig(NeXT.Matrices.A));
        imag_part_NeXT=imag(eig(NeXT.Matrices.A));

        
        %Plot OMA
        %set(gcf,'units','points','position',[100,100,700,700])
        subplot(2,2,1)
        circle(0,0,1)
        hold on

        for i=1:opt_order
            quiver (0,0,real_part_OMA(i),imag_part_OMA(i), 'r')
            hold all
        end
        xlim([-1.25 1.25])
        ylim([-1.25 1.25])
        grid on
        title(['OMA-Covar for Subject: ' num2str(subject)])
        %filename=[basename, num2str(trial) ,extension]

        % OMA-data plot
        subplot(2,2,2)
        circle(0,0,1)
        hold on

        for i=1:opt_order
            quiver (0,0,real_part_data(i),imag_part_data(i), 'b')
            hold all
        end
        xlim([-1.25 1.25])
        ylim([-1.25 1.25])
        grid on
        title(['OMA-Data for Subject: ' num2str(subject)])


        % Plot n4sid
        subplot(2,2,3)
        circle(0,0,1)
        hold on

        for i=1:length(real_part_n4)
            quiver (0,0,real_part_n4(i),imag_part_n4(i), 'k')
            hold all
        end
        xlim([-1.25 1.25])
        ylim([-1.25 1.25])
        grid on
        title(['n4sid Subject: ' num2str(subject)])
        
        %Plot NeXT-ERA
        subplot(2,2,4)
        circle(0,0,1)
        hold on

        for i=1:length(real_part_NeXT)
            quiver (0,0,real_part_NeXT(i),imag_part_NeXT(i), 'Color', loyolagreen)
            hold all
        end
        xlim([-1.25 1.25])
        ylim([-1.25 1.25])
        grid on
        title(['ERA Subject: ' num2str(subject)])
        
        %close all
    end
    sgtitle('Pole Comparison of Output Only SID Algorithms for Subject Over All Trials')

    saveas(gcf,filename)
    clear s02
    close all
end