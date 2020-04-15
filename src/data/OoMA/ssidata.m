function [A,C,G,R0] = ssidata(Y,order,s)
% [A,C,G,R0] = SSIDATA(Y,order,s)
%   Data-based stochastic subspace identification (SSI-data).
%
%   INPUTS:
%	Y       sensor data matrix
%   order	desired maximum model order to identify (scalar)
%   s       number of block rows in the block hankel matrix, should be at 
%           least ceil(order/ns) to obtain results up to the desired model
%           order, generally recommended that s > order
%
%   OUTPUTS:
%   A       cell array of state transition matrices for model orders {i}
%   C       cell array of output matrices for model orders {i}
%   G       cell array of next state output covariance matrices for model
%           orders {i}
%   R0      zero-lag output covariances
%
% REFERENCES:
% [1]   "Subspace Identification for Linear Systems" by Peter van Overschee
%       and Bart de Moor, doi:10.1007/978-1-4613-0465-4
% [2]	"System Identification Methods for (Operational) Modal Analysis:
%       Review and Comparison" by Edwin Reynders in Archives of
%       Computational Methods in Engineering, Vol. 19, No. 1,
%       doi:10.1007/s11831-012-9069-x
% [3]   "Operational Modal Analysis of Civil Engineering Structures" by
%       Carlo Rainieri and Giovanni Fabbrocino,
%       doi:10.1007/978-1-4939-0767-0
%
% NOTES:
% (1)   this implementation can easily consume a lot of memory for large
%       sensor records or numerous time lags, iterative algorithms
%       operating on the sensor data to create the projection are
%       recommended and will be the focus of future updates

disp('SSI-data status:')

[r,c] = size(Y);
if r > c % make sure y is shaped correctly with samples going across rows
    Y = Y';
end

[ns,nt] = size(Y); % ns = # of sensors, nt = # of samples

% shifted data matrix
disp('  forming shifted data matrix...')
Yh = zeros(ns*2*s,nt-2*s+1);
for i = 1:2*s % go down block rows of the Hankel data matrix
    Yh((i-1)*ns+1:i*ns,:) = Y(:,i:nt-2*s+i); % fill out the entire row
end
Yh = Yh/sqrt(nt-2*s+1);

% QR decomposition and projection of raw data
disp('  projecting raw data...')
R = triu(qr(Yh'))';
R = R(1:2*s*ns,1:2*s*ns);
Proj = R(ns*s+1:2*ns*s,1:ns*s);

% SVD (no weighting = balanced PC)
disp('  performing singular value decomposition...')
[U,S,~] = svd(Proj);
S = diag(S);

% zero lag output covariance
R0 = R(ns*s+1:ns*(s+1),:)*R(ns*s+1:ns*(s+1),:)'; 

% output cell arrays
A = cell(1,order);
C = cell(1,order);
G = cell(1,order);

% loop over model orders and generate system matrices
disp(['  generating system matrices A,C,G for ' num2str(order) ' model orders...'])
for i = 1:order
    U1 = U(:,1:i);
    gam  = U1*diag(sqrt(S(1:i)));
    gamm = U1(1:ns*(s-1),:)*diag(sqrt(S(1:i)));
    gam_inv  = pinv(gam);
    gamm_inv = pinv(gamm);
    A{i} = gamm_inv*gam(ns+1:ns*s,:); % state transition matrix
    C{i} = gam(1:ns,:); % output matrix
    delta = gam_inv*(R(ns*s+1:2*ns*s,1:ns*s)*R(1:ns*s,1:ns*s)');
    G{i} = delta(:,ns*(s-1)+1:ns*s); % next state output covariance matrix
end

disp('SSI-data finished.')


end