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

%% n4sid
data = iddata(Y',[],T);
sys = n4sid(data,10);
[A_n4,~,C_n4,~] = ssdata(sys);

