function mm = macmatrix(Phi1,Phi2)
% mm = MACMATRIX(Phi1,Phi2)
%   Returns a matrix containing the modal assurance criterion (MAC) values
%   between two mode shape matrices. Mode shape vectors are oriented as
%   column vectors within the input matrices.
%
%   INPUTS:
%   Phi1    first matrix of mode shape vectors
%   Phi2    second matrix of mode shape vectors
%
%   OUTPUTS:
%   mm      a matrix containing the modal assurance criterion values
%           between the two matrices of mode shape vector where mm(i,j)
%           is the MAC value between Phi1(:,i) and Phi2(:,j)

% assumes mode shape vectors are row vectors
nm_1 = size(Phi1,2);
nm_2 = size(Phi2,2);

mm = zeros(nm_1,nm_2);

for i = 1:nm_1
    for j = 1:nm_2
        mm(i,j) = mac(Phi1(:,i),Phi2(:,j));
    end
end


end