[r,c] = size(Y_final);
if r > c % make sure y is shaped correctly with samples going across rows
    Y_final = Y_final';
end

[ns,nt] = size(Y_final)
s=6;
order=25;
disp('  forming shifted data matrix...')
Yh = zeros(ns*2*s,nt-2*s+1)
for i = 1:2*s % go down block rows of the Hankel data matrix
    Yh((i-1)*ns+1:i*ns,:) = Y_final(:,i:nt-2*s+i); % fill out the entire row
end
Yh = Yh/sqrt(nt-2*s+1)

disp('  projecting raw data...')
R = triu(qr(Yh'))';
R = R(1:2*s*ns,1:2*s*ns);
Proj = R(ns*s+1:2*ns*s,1:ns*s);