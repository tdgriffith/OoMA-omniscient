%% Generate an A matrix, compare timeseries results to original data
s02=load('s09.mat');
subject=09;
basename='S015E';
extension='.png';
trial=1;

Y=s02.data(trial,1:6,:);
Y=squeeze(Y);
t = linspace(0,63,length(Y));
fs=128;
T=1/fs;

[Syy,freqs] = pwelch(Y',[],[],[],fs); % obtain estimates of the output power spectrums
clf
subplot(2,1,1)
plot(t,Y)
xlabel('Time (s)')
ylabel('Acceleration (ms^{-2})')
axis tight
subplot(2,1,2)
plot(freqs,10*log10(Syy))
xlabel('Frequency (Hz)')
ylabel('PSD (dB)')
axis tight