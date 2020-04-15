function [A,C] = ssicovref(Y,order,s)
% [A,C] = SSICOVREF(Y,order,S)
%	Covariance-based stochastic subspace identification with multiple
%   setups and reference sensors (SSI-cov/ref).
%
%   INPUTS:
%   Y       structure array containing the reference (Y(i).ref) and moving
%           (Y(i).mov) sensor data for sensor setup i, time records must be
%           equal in length but the number of moving sensors can vary
%           between sensor setups
%   order	maximum model order for the identified system, scalar
%   s    	number of time lags used in formation of the block Hankel
%           matrix, as a rule of thumb s > order
%
%   OUTPUTS:
%   A     	cell array of state transition matrices for model orders {i}
%   C     	cell array of output matrices for model orders {i}
%
% REFERENCES:
% [1]   "Modular Subspace-Based System Identification from Multi-Setup
%       Measurements" by Michael Dohler and Laurent Mevel in IEEE
%       Transactions on Automatic Control, Vol. 57, No. 11,
%       doi:10.1109/TAC.2012.2193711
%
% NOTES:
% (1)	memory requirements could be reduced by iteratively forming the [A]
%       and [C] matrices as suggested by the referenced work
% (2)   the setup used for state-basis scaling is simply chosen as the
%       first sensor setup; alternately, condition number checking could be
%       used as a selection criteria
% (3)   search for a method to form the next-state output covariance matrix
%       [G] and zero lag output covariances [R0] for ssi-cov/ref so output
%       spectrums can be generated as well as solution of a forward
%       innovation model for perfoming model synthesis from selected system
%       poles
% (4)   no input checking is performed, will be added in future updates

disp('SSI-cov/ref status:')

n_setups = numel(Y); % number of different sensor arrangements/setups
n_ref = size(Y(1).ref,1); % number of reference sensors

% number of moving sensors in each setup
n_mov = zeros(1,n_setups); 
for i = 1:n_setups
    n_mov(i) = size(Y(i).mov,1);
end

n_s = n_ref+sum(n_mov); % total number of sensors

Obs = cell(n_setups,1);

for i = 1:n_setups
    disp(['  processing setup ' num2str(i) ' of ' num2str(n_setups) '...'])
    % generate hankel matrix estimate
    H = blockhankel(Y(i).ref,[Y(i).ref; Y(i).mov],s);
    % SVD the hankel matrix
    [U,S,~] = svd(H,0);
    % truncate at maximum model order
    S = S(1:order,1:order);
    % obtain the observability matrix for this setup
    obs = U(:,1:order)*sqrt(S);
    % indexing vectors to get block rows without looping
    id = zeros(n_ref,s);
    for j = 1:n_ref
        id(j,:) = (0:s-1)*(n_ref+n_mov(i))+j;
    end
    % obtain the reference portion of the observability matrix
    obs_ref = obs(id,:);
    obs(id,:) = []; % delete reference portion to from obs_mov
    if i == 1
        % state basis for scaling, just using the first setup
        % could check for conditioning to choose which setup
        obs1_ref = obs_ref;
    end
    % scale the obs matrix to the global state basis
    obs = obs*pinv(obs_ref)*obs1_ref;
    Obs{i} = obs;
end

% global observation matrix formation via block-interleaving
disp('  generating global observability matrix...')
Obs_all = zeros(n_s*s,order);
for i = 1:s
    % reference portion block rows
    id1 = (i-1)*n_s+1;
    id2 = id1+n_ref-1;
    Obs_all(id1:id2,:) = obs1_ref((i-1)*n_ref+1:i*n_ref,:);
    for j = 1:n_setups
        % moving sensor portion block rows
        id1 = id2+1;
        id2 = id1+n_mov(j)-1;
        Obs_all(id1:id2,:) = Obs{j}((i-1)*n_mov(j)+1:i*n_mov(j),:);
    end
end
clear Obs

% output cell arrays
A = cell(order,1);
C = cell(order,1);

% loop over model orders
disp(['  generating system matrices A,C for ' num2str(order) ' model orders...'])
for i = 1:order
    A{i} = pinv(Obs_all(1:end-n_s,1:i))*(Obs_all(n_s+1:end,1:i)); % system matrix
    C{i} = Obs_all(1:n_s,1:i); % output matrix
end

disp('SSI-cov/ref finished.')


end