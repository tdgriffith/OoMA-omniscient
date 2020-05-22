% load(['ecog_window.mat']);
mass=[4,0;0,3];
stiff=[2,-1;-1,1];
damp=[0.4,-0.1; -0.1, 0.1];
A=[zeros(length(mass)),eye(length(mass));-inv(mass)*stiff,-inv(mass)*damp];
[eigV,eigD]=eig(A);

sys=ss(A,[1;1;0;0],eye(4),[0;0;0;0]);
[X,t]=step(sys);
X=X';
X=detrend(X);
% parameters:
r = 10; % number of modes
nstacks = 5; % number of stacks
data_org=X;
% construct the augmented, shift-stacked data matrices
Xaug = [];
for st = 1:nstacks,
    Xaug = [Xaug; X(:, st:end-nstacks+st)];%#ok<AGROW>
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
dt=t(2)-t(1)
omega = diag(log(diag(lambda)))/dt
Omega=diag(omega);
%% eigenvalue
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
axis([-2 2 -2 +2]);
axis square;

set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6 3], 'PaperPositionMode', 'manual');
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
xlim([0 150]);
axis square;

% power spectrum
timesteps = size(X, 2);
srate = 1/dt;
nelectrodes = 59;
NFFT = 2^nextpow2(timesteps);
f = srate/2*linspace(0, 1, NFFT/2+1);

subplot(1,2,2); 
hold on;
for c = 1:nelectrodes,
    fftp(c,:) = fft(X(c,:), NFFT);
    plot(f, 2*abs(fftp(c,1:NFFT/2+1)), ...
        'Color', 0.6*[1 1 1]);
end;
plot(f, 2*abs(mean(fftp(c,1:NFFT/2+1), 1)), ...
    'k', 'LineWidth', 2);
xlim([0 150]);
ylim([0 400]);
axis square;
box on;

set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6 3], 'PaperPositionMode', 'manual');
% print('-depsc2', '-loose', ['../figures/DMD_spectrum.eps']);

c=7;
Phi_phys=Phi(1:c,:);
Phi_max=max(Phi_phys);
Phi_phys_norm=Phi_phys./Phi_max;

figure
plot(t(1:150),data_org(1,1:150))
%%
b = Phi\Xaug(:,1);
for k=1:length(t)
    time_dyanmics(:,k)=(b.*exp(omega*t(k)));
end
Xdmd=Phi*time_dyanmics;