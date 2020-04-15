%% Setup
for subject = 1:32
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
    basename=['S' num2str(subject1) 'S' num2str(subject2) '-compE'];
    extension='.png';
    s02=load(load_name2)

%% Loop Over Trials
    for trial = 1:40
        % Trial Data
        Y1=s01.data(trial,1:6,:);
        Y1=squeeze(Y1);
        Y2=s02.data(trial,1:6,:);
        Y2=squeeze(Y2);
        t = linspace(0,63,length(Y1));
        fs=128;
        T=1/fs;
        % TS Plots, don't need right now
        %[Syy,freqs] = pwelch(Y1',[],[],[],fs); % obtain estimates of the output power spectrums
        %[Syy2,freqs2] = pwelch(Y2',[],[],[],fs); % obtain estimates of the output power spectrums

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

        % OMA Covariance Algorithm
        order = 25;
        s = 2*order;
        opt_order=12;
        [A_cov,C_cov,G_cov,R0_cov] = ssicov(Y1,order,s);
        [A_cov2,C_cov2,G_cov2,R0_cov2] = ssicov(Y2,order,s); 
        eig(A_cov{opt_order})
        err = [0.01,0.05,0.98];
        %[IDs_cov] = plotstab(A_cov,C_cov,Y,T,[],err);
        [fn_cov,zeta_cov,Phi_cov] = modalparams(A_cov,C_cov,T);
        [fn_cov2,zeta_cov2,Phi_cov2] = modalparams(A_cov2,C_cov2,T);
        real_part=real(eig(A_cov{opt_order}));
        imag_part=imag(eig(A_cov{opt_order}));
        real_part2=real(eig(A_cov2{opt_order}));
        imag_part2=imag(eig(A_cov2{opt_order}));

        %n4sid
        data1 = iddata(Y1',[],T);
        data2 = iddata(Y2',[],T);
        sys1 = n4sid(data1,12);
        sys2 = n4sid(data2,12);
        [A_n4_1,~,C_n4_1,~] = ssdata(sys1);
        [A_n4_2,~,C_n4_2,~] = ssdata(sys2);
        real_part_n4=real(eig(A_n4_1));
        imag_part_n4=imag(eig(A_n4_1));
        real_part_n4_2=real(eig(A_n4_2));
        imag_part_n4_2=imag(eig(A_n4_2));

        % OMA-data
        [A_data,C_data,G_data,R0_data] = ssidata(Y1,order,s);
        real_part_data=real(eig(A_data{opt_order}));
        imag_part_data=imag(eig(A_data{opt_order}));
        [A_data2,C_data2,G_data2,R0_data2] = ssidata(Y2,order,s);
        real_part_data2=real(eig(A_data2{opt_order}));
        imag_part_data2=imag(eig(A_data2{opt_order}));

        %NeXT-ERA
        [NeXT]=NExTFERA(Y1,6,2000,4,0.1,fs,800,200,12,10,1);
        real_part_NeXT=real(eig(NeXT.Matrices.A));
        imag_part_NeXT=imag(eig(NeXT.Matrices.A));
        [NeXT2]=NExTFERA(Y2,6,2000,4,0.1,fs,800,200,12,10,1);
        real_part_NeXT2=real(eig(NeXT2.Matrices.A));
        imag_part_NeXT2=imag(eig(NeXT2.Matrices.A));

        figure
        set(gcf,'units','points','position',[500,-300,700,700])

        subplot(2,2,1)
        circle(0,0,1)
        hold on

        for i=1:8
            quiver (0,0,real_part(i),imag_part(i), 'r')
            hold all
            quiver (0,0,real_part2(i),imag_part2(i), 'b')        
            hold all
        end
        title('OMA-Covar')
        xlim([-1.25 1.25])
        ylim([-1.25 1.25])
        grid on

        % OMA-data plot
        subplot(2,2,2)
        circle(0,0,1)
        hold on

        for i=1:8
            quiver (0,0,real_part_data(i),imag_part_data(i), 'r')
            hold all
            quiver (0,0,real_part_data2(i),imag_part_data2(i), 'b')
            hold all
        end
        xlim([-1.25 1.25])
        ylim([-1.25 1.25])
        grid on
        title('OMA-Data')

        % Plot n4sid
        subplot(2,2,3)
        circle(0,0,1)
        hold on

        for i=1:8
            quiver (0,0,real_part_n4(i),imag_part_n4(i), 'r')
            hold all
            quiver (0,0,real_part_n4_2(i),imag_part_n4_2(i), 'b')
            hold all
        end
        xlim([-1.25 1.25])
        ylim([-1.25 1.25])
        grid on
        title('n4sid')

        %Plot NeXT-ERA
        subplot(2,2,4)
        circle(0,0,1)
        hold on

        for i=1:8
            quiver (0,0,real_part_NeXT(i),imag_part_NeXT(i), 'r')
            hold all
            quiver (0,0,real_part_NeXT2(i),imag_part_NeXT2(i), 'b')
            hold all
        end
        xlim([-1.25 1.25])
        ylim([-1.25 1.25])
        grid on
        title(['ERA Subject: ' num2str(subject)])


        sgtitle(join((['Subject: ' num2str(subject1), num2str(subject2) ' Emotion: ' onlineratings(trial)])))
        filename=[basename, num2str(trial) ,extension]
        saveas(gcf,filename)
        close all
    end
end