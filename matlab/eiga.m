function [U,D]=eiga(A)
% EIGA normalized eigenvalue decomposition
%     It is basically equivalent to eig, but the eigenvectors are sorted 
%     in terms of the eigenvalues (in accending order), and have 
%     nonnegative first entries. 

[U,D]=eig(A);
dd=diag(D);
[~,isD]=sort(dd);
D=diag(dd(isD));
sgnflip=sign(U(1,:)); sgnflip(sgnflip==0)=1;
U=U(:,isD); U=U*diag(sgnflip);