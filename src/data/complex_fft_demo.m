Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 6000;             % Length of signal
t = (0:L-1)*T;        % Time vector
S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
X = S + 2*randn(size(t));
plot(1000*t(1:50),X(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P3=Y(1:L/2+1);
P3norm=P3./norm(P3);
P4=abs(P3);
Phase4=angle(P3);
f = Fs*(0:(L/2))/L;
figure
subplot(1,2,1)
plot(f,P4./max(abs(P4)))
subplot(1,2,2)
plot(f,P1./max(abs(P1)))


figure
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
%% Compare with small window: Nyquist and whatnot
X2=X(:,1:L/8);
L2=size(X2,2);
Y2=fft(X2);
Pshort=abs(Y2/L2);
Pshort2 = Pshort(1:L2/2+1);
Pshort2(2:end-1) = 2*Pshort2(2:end-1);
f2 = Fs*(0:(L2/2))/(L2);
figure
plot(f2,Pshort2(1,:)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% Reduce order: Conclusion, n on fft doesn't reduce basis space, limits data points
% X2=X(:,1:L/8);
% L2=size(X2,2);
% Y2=fft(X2,150);
% L2=size(Y2,2);
% Pshort=abs(Y2/L2);
% Pshort2 = Pshort(1:L2/2+1);
% Pshort2(2:end-1) = 2*Pshort2(2:end-1);
% f2 = Fs*(0:(L2/2))/(L2);
% figure
% plot(f2,Pshort2(1,:)) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
%% DCT
I = imread('cameraman.tif');
I = im2double(I);

T = dctmtx(8);
dct = @(block_struct) T * block_struct.data * T';
B = blockproc(I,[8 8],dct);

mask = [1   1   1   1   0   0   0   0
        1   1   1   0   0   0   0   0
        1   1   0   0   0   0   0   0
        1   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0];
B2 = blockproc(B,[8 8],@(block_struct) mask .* block_struct.data);

invdct = @(block_struct) T' * block_struct.data * T;
I2 = blockproc(B2,[8 8],invdct);

imshow(I)
figure
imshow(I2)

%% Image but with DFT
I = imread('cameraman.tif');
I = im2double(I);

T = dftmtx(8);
dct = @(block_struct) T * block_struct.data * T';
B = blockproc(I,[8 8],dct);

mask = [1   1   1   1   0   0   0   1
        1   1   1   0   0   0   0   0
        1   1   0   0   0   0   0   0
        1   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        1   0   0   0   0   0   0   1
        1   1   0   0   0   0   1   1];
%mask=eye(8);
B2 = blockproc(B,[8 8],@(block_struct) mask .* block_struct.data);

invdct = @(block_struct) T' * block_struct.data * T;
I2 = blockproc(B2,[8 8],invdct);
I2=real(I2);
I2=I2./(length(T)^2);
imshow(I)
figure
imshow(I2)

%% Reduce order: Use DFT matrix
X2=X(:,1:L/4);
%Y2= dftmtx(size(X2,2))*X2';
Y2= dftmtx(40)*X2';
Y2=Y2';
L2=size(X2,2);

Pshort=abs(Y2/L2);
Pshort2 = Pshort(1:L2/2+1);
Pshort2(2:end-1) = 2*Pshort2(2:end-1);
f2 = Fs*(0:(L2/2))/(L2);
figure
plot(f2,Pshort2(1,:)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')