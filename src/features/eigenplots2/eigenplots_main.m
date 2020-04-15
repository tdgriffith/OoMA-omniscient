
s02=load('s09.mat');
subject=09;
basename='S015E';
extension='.png';
figure
for trial = 1:40

    Y=s02.data(trial,1:6,:);
    Y=squeeze(Y);
    t = linspace(0,63,length(Y));
    fs=128;
    T=1/fs;

    [Syy,freqs] = pwelch(Y',[],[],[],fs); % obtain estimates of the output power spectrums
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
    order = 25;
    s = 2*order;
    opt_order=12;
    [A_cov,C_cov,G_cov,R0_cov] = ssicov(Y,order,s);
    eig(A_cov{opt_order})
    err = [0.01,0.05,0.98];
    %[IDs_cov] = plotstab(A_cov,C_cov,Y,T,[],err);
    [fn_cov,zeta_cov,Phi_cov] = modalparams(A_cov,C_cov,T);

    real_part_OMA=real(eig(A_cov{opt_order}));
    imag_part_OMA=imag(eig(A_cov{opt_order}));
    
    %n4sid
    data = iddata(Y',[],T);
    sys = n4sid(data,12);
    [A_n4,~,C_n4,~] = ssdata(sys);
    real_part=real(eig(A_n4));
    imag_part=imag(eig(A_n4));    

    %figure
    set(gcf,'units','points','position',[200,200,400,400])
    circle(0,0,1)
    hold on

    for i=1:8
        quiver (0,0,real_part(i),imag_part(i), 'r')
        hold all
    end
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])
    grid on
    title(['Subject: ' num2str(subject)])
    filename=[basename, num2str(trial) ,extension]
    %saveas(gcf,filename)
    %close all
end