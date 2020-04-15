function [SV,F,U] = cmif(Y,win,dt)
% [SV,F] = CMIF(Y,win,dt)
%   Complex Mode Indicator Function (CMIF). Returns the singular values of
%   the cross power spectrums as a function of frequency.
%
%   INPUTS:
%   Y       matrix containing output samples
%   win     (optional) window to use for estimation of the cross PSD
%   dt      sampling period of the measured system
%
%   OUTPUTS:
%   SV      matrix containing the singular values of the cross PSD as a
%           function of frequency where SV(n,:) is the nth singular value
%           for all frequencies up to the nyquist frequency (1/dt/2)
%   F       frequency vector (in Hz) that corresponds to the second
%           dimension of SV
%
% REFERENCES:
% [1]   Brincker, Rune, Lingmi Zhang, and P. Andersen. "Modal
%       identification from ambient responses using frequency domain
%       decomposition." Proc. of the 18th International Modal Analysis
%       Conference (IMAC), San Antonio, Texas. 2000.
    
% input conditioning, samples go down rows (assumes more samples than
% channels)
[r,c] = size(Y);
if r < c
    Y = Y';
end

if isempty(dt)
    dt = 1;
end

[~,ns] = size(Y); % number of output states

% calculate the cross PSD using CPSD
[Syy,F] = cpsd(Y,Y,win,[],[],1/dt,'mimo');

SV = zeros(ns,length(F));
U = zeros(ns,ns,length(F));
% perform SVD at each frequency point
for i = 1:length(F)
    [U(:,:,i),S,~] = svd(squeeze(Syy(i,:,:)));
    SV(:,i) = diag(S);
end


end