function [f,zeta,Phi] = modalparams(A,C,dt)
% [f,psi,Phi] = MODALPARAMS(A,C,dt)
%	Modal decomposition of discrete state space system.
%
%   INPUTS:
%   A       cell array of system matrices for model order {i}
%   C       cell array of output matrices for model order {i}
%   dt      sampling period of the discrete system A,C
%
%   OUTPUTS:
%   f       cell array containing the system pole frequencies in Hz for
%           model orders {i}
%   psi     cell array containing the damping ratios of each pole for model
%           orders {i}
%   Phi     cell array containing the mode shape vectors for model orders
%           {i}

% NOTES:
% (1)	modal scaling (normalization) is not performed
% (2)   complex conjugate pairs are eliminated and modes are sorted by
%       frequency

% check for cell array input
if ~iscell(A)
    A = {A};
end
if ~iscell(C)
    C = {C};
end

f = cell(size(A));
zeta = cell(size(A));
Phi = cell(size(A));

% loop over model orders
for i = 1:length(A)
    [v,d] = eig(A{i});
    lam = log(diag(d))/dt;
    f{i} = abs(lam)/2/pi; % modal frequencies (Hz)
    [f{i},I] = sort(f{i}); % sort using ascending frequencies
    zeta{i} = -real(lam)./abs(lam); % modal damping ratios
    zeta{i} = zeta{i}(I);
    Phi{i} = C{i}*v; % mode shapes
    Phi{i} = Phi{i}(:,I);
    % eliminate complex conjugate pairs
    [f{i},I] = unique(f{i});
    zeta{i} = zeta{i}(I);
    Phi{i} = Phi{i}(:,I);
end



end