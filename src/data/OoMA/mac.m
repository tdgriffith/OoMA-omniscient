function m = mac(phi1,phi2)
% m = MAC(phi1,phi2)
%	Returns the modal assurance criterion (MAC) between two mode shape
%	vectors.
%
%	INPUTS:
%	phi1	first mode shape vector
%   phi2	second mode shape vector
%
%   OUTPUTS:
%   m       the modal assurance criterion for the two mode shape vectors
%
% REFERENCES:
% [1]   Allemang, Randall J. "The modal assurance criterion–twenty years of
%       use and abuse." Sound and vibration 37.8 (2003): 14-23.

m = (abs(phi1'*phi2))^2/((phi1'*phi1)*(phi2'*phi2));


end