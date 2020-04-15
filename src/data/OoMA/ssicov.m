function [A,C,G,R0] = ssicov(Y,order,s)
% [A,C,G,R0] = SSICOV(Y,order,s)
%	Covariance-based stochastic subspace identification (SSI-cov).
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
% [1]	"System Identification Methods for (Operational) Modal Analysis:
%       Review and Comparison" by Edwin Reynders in Archives of
%       Computational Methods in Engineering, Vol. 19, No. 1,
%       doi:10.1007/s11831-012-9069-x
% [2]   "Operational Modal Analysis of Civil Engineering Structures" by
%       Carlo Rainieri and Giovanni Fabbrocino,
%       doi:10.1007/978-1-4939-0767-0

disp('SSI-cov status:')

% input conditioning, assumes more samples than sensors
[r,c] = size(Y);
if r > c % make sure y is shaped correctly with samples going across rows
    Y = Y';
end
[ns,nt] = size(Y); % ns = # of sensors, nt = # of samples

% set the number of block rows in the block hankel matrix
if isempty(s) || s < ceil(order/ns)
    s = ceil(order/ns); % s is at least ceil(order/ns), better results obtained using more time lags
    disp(['warning, block hankel matrix size too small, using s = ' num2str(ceil(order/ns))])
end

R0 = 1/(nt)*Y(:,1:nt)*Y(:,1:nt)'; % zero lag output covariances

% block Hankel matrix of output covariances
disp('  forming block Hankel matrix...')
H = blockhankel(Y,Y,s);

% balanced realization (no weighting matrices)
disp('  performing singular value decomposition...')
[U,S,V] = svd(H,0); % decomposition

% truncate at max model order
S = S(1:order,1:order);
S = diag(S);
U = U(:,1:order);
V = V(:,1:order);

% output cell arrays
A = cell(order,1);
C = cell(order,1);
G = cell(order,1);

% loop over model orders
disp(['  generating system matrices A,C,G for ' num2str(order) ' model orders...'])
for i = 1:order
    U1 = U(:,1:i); % truncate decomposition at model order
    V1 = V(:,1:i); % truncate decomposition at model order
    ss = diag(sqrt(S(1:i))); % square root of truncated singular values
    Obs = U1*ss; % observability matrix
    Con = ss*V1'; % controllabiliy matrix (reversed if using Toeplitz)
    A{i} = pinv(Obs(1:end-ns,:))*(Obs(ns+1:end,:)); % system matrix
    C{i} = Obs(1:ns,:); % output matrix
    G{i} = Con(:,1:ns); % next state output covariance matrix
    % G{i} = Con(:,end-ns+1:end); % if using Toeplitz matrix
end

disp('SSI-cov finished.')


end