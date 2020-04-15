function [Syy,F] = ssispectrum(A,C,G,R0,dt,order,n)
% [Syy,F] = SSISPECTRUM(A,C,G,R0,dt,n)
%   Stochastic system power spectrum (normalized).
%
%   INPUTS:
%   A       discrete state transition matrix
%   C       output matrix
%   G       next-state output matrix
%   R0      zero lag output covariance matrix
%   dt      sampling period of the discrete system
%   n       number of points to use in formation of the output power
%           spectrum, up to the nyquist frequency
%
%   OUTPUTS:
%   Syy     3-dimensional array containing the output power spectrums,
%           where Syy(i,j,:) is the power spectrum between outputs i and j
%   F       frequency values for the output power spectrums, in Hz
%
% REFERENCES:
% [1]   "Operational Modal Analysis of Civil Engineering Structures" by
%       Carlo Rainieri and Giovanni Fabbrocino,
%       doi:10.1007/978-1-4939-0767-0
%
% NOTES:
% (1)   output normalization is performed using the largest magnitude
%       component

fs = 1/dt;
[no,ns] = size(C{order});
I = eye(ns);
F = linspace(0,fs/2,n);
z = exp(1i*F*2*pi*dt);
Syy = zeros(no,no,length(z));

for i = 1:length(z)
    Syy(:,:,i) = C{order}/(z(i)*I-A{order})*G{order}+R0+G{order}'/(1/z(i)*I-A{order}')*C{order}';
end

Syy = Syy/max(abs(Syy(:)));


end