function IDs = plotstab(A,C,Y,dt,win,err)
% IDs = PLOTSTAB(f,psi,y,dt,err)
%   Creates a stabilzation diagram for modal analysis purposes in the
%   current figure.
%
%   INPUTS:
%   A       cell array of system matrices
%	C       cell array of output matrices
%	Y       test data used for model identification purposes
%	dt      sampling period of output data y
%	win     optional (replace with []) window to be used for estimation of
%           output power spectrums
%	err     3-element vector of percent errors for stability criteria
%           (frequency, damping, and modal assurance criterion),
%           default = [0.01,0.05,0.98]
%
%	OUTPUTS:
%	IDs     cell array containing logical vectors of indices of stable
%           modes on diagram for all model orders of the identified system

% check if user specified err vector
if isempty(err)
    warning('No stabilization criteria specified, using default settings for stabilization criteria')
    err = [0.01 0.05 0.98];
end

% check for out of bounds on error vector
if any(err<0) && any(err>1)
    error('Vector err out of bounds in the range 0 to 1')
end



% input conditioning, samples go down rows (assumes more samples than
% channels)
[r,c] = size(Y);
if r < c
    Y = Y';
end

% generate the complex mode indicator function (CMIF)
[SV,F] = cmif(Y,win,dt);

% plot the first singular value of the CMIF to aid in pole selection
clf
yyaxis('right')
plot(F,10*log10(SV(1,:)))

% generate modal decompositions
[f,psi,Phi] = modalparams(A,C,dt);

% loop over model orders
IDs = cell(size(A));
for i = 1:length(A)-1
    [f1,I1,~] = unique(f{i}); % ignore repeated modes
    [f2,I2,~] = unique(f{i+1}); % ignore repeated modes
    psi1 = psi{i}(I1);
    psi2 = psi{i+1}(I2);
    phi1 = Phi{i}(:,I1);
    phi2 = Phi{i+1}(:,I2);
    % frequency stability criteria
    ef = logical(sum(abs(bsxfun(@minus,f1,f2')./f1)<=err(1),2));
    % damping stability criteria
    epsi = logical(sum(abs(bsxfun(@minus,psi1,psi2')./psi1)<=err(2),2));
    % MAC stability criteria
    mac_vals = zeros(1,length(f2));
    ephi = zeros(length(f1),1);
    % check each mode shape vector with ever other mode shape vector from
    % the next model order up
    for j = 1:length(f1)
        for k = 1:length(f2)
            mac_vals(k) = mac(phi1(:,j),phi2(:,k));
        end
        ephi(j) = logical(sum(mac_vals>=err(3)));
    end
    % valid (stable) poles
    IDs{i} = I1(ef&epsi&ephi);
    stable = f1(ef&epsi&ephi);
    unstable = f1(~(ef&epsi&ephi));
    % add unstable poles to plot
    if ~isempty(unstable)
        yyaxis('left')
        a = plot(unstable,i,'.');
        hold on
    end
    % add stable poles to plot
    if ~isempty(stable)
        yyaxis('left')
        hold on
        b = plot(stable,i,'x');
    end
end

% label axes
title('Stabilization Diagram')
xlabel('Frequency (Hz)')
yyaxis('left')
ylabel('Model Order')
yyaxis('right')
ylabel('CMIF Magnitude (dB)')
xlim([F(1) F(end)])
legend([a(1); b(1)],'unstable','stable')
hold off


end