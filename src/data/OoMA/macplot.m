function macplot(mm)
% MACPLOT(mm)
%   Plots the MAC matrix mm in the current figure.
%   
%   INPUTS:
%   mm      the matrix containing the MAC values obtained from calling
%           macmatrix on two sets of mode shape vectors. Mode shape vectors
%           are column vectors.

clf
imagesc(mm')
axis xy
colorbar
title('Modal Assurance Criterion')
xticks(1:size(mm,1))
yticks(1:size(mm,2))


end