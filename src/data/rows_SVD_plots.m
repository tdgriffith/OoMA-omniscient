%% Singular value plots for OMA
disp('SSI-data status:')
load_name1=['S02_restingPre_EO.mat']
s01=load(load_name1);
Y1=s01.dataRest(1:64,:);
Y=squeeze(Y1);
order=40;
s_out={};
count=0;
for s=10:100:1310
    count=count+1;

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
    s_out{count}=S;
    
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
    count/38
end
cutoff=length(s_out{1});
for i=1:size(s_out,2)
    s_out{i}=s_out{i}(1:cutoff);
end
s_out_arr=cell2mat(s_out);




