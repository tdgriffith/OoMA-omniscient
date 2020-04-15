function H = blockhankel(Y_ref,Y_all,s)
% H = BLOCKHANKEL(Y_ref,Y_all,s)
%	Generate block Hankel matrix using covariances from output data.
%
%	INPUTS:
%	Y_ref   reference data
%	Y_all   moving sensor data, if not using reference-based
%           identification then Y_all = Y_ref
%   s     	number of time lags used in covariance calculation
%
%   OUTPUTS:
%   H       block Hankel matrix

r = size(Y_all,1);
c = size(Y_ref,1);
R = zeros(r,c,2*s);
% loop over sensor channels
for i = 1:r
    for j = 1:c
        % use FFT-based, built-in xcov() to estimate output covariances
        temp = xcov(Y_all(i,:),Y_ref(j,:),2*s,'biased');
        R(i,j,:) = temp(2*s+2:end); % only the positive time lags
    end
end

% form the block Hankel matrix in rows
for i = 1:s
    temp = reshape(R(:,:,i:s+i-1),r,c*s); % one block row at time lag i
    H((i-1)*r+1:i*r,:) = temp; % Hankel matrix
    % T(:,ns*s-i*ns+1:ns*s-(i-1)*ns) = temp'; % Toeplitz matrix (for reference)
end


end