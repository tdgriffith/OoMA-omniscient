function [Pxx,F,tvec]=compute_eeg_power(x,Fs,tf)
    dt=1/Fs;
    t = 0:dt:tf-dt;
    winsize = 150;
    numoverlap = round(0.85*winsize);
    win = hamming(winsize);
    X = buffer(x,winsize,numoverlap);
    for nn = 1:size(X,2)
    [Pxx(:,nn),F] = pwelch(X(:,nn),win,length(win)/2,length(win),Fs);
    end
     % create a time vector
     idxbegin = find(X(:,1) == 0);
     numpresteps = length(idxbegin);
     idxend = find(X(:,end) == 0);
     numpoststeps = length(idxend);
     tbegin = -(numpresteps*dt)/2;
     tend = t(end)+((numpoststeps*dt))/2;
     tvec = linspace(tbegin,tend,size(Pxx,2));
%      surf(tvec,F,10*log10(abs(Pxx)),'EdgeColor','none');   
%      axis xy; axis tight; colormap(jet); view(0,90);
%      xlabel('Time (sec)');
%      ylabel('Frequency (Hz)');
end