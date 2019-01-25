function [W,L,U]=skeleton15(isnorm,ws)
if nargin<1
    isnorm=0;
end
alist=[3,2,2,2,4,7,5,8,3,3,10,13,11,14;
       2,1,4,7,5,8,6,9,10,13,11,14,12,15];
if nargin<2
    wlist=ones(1,14);
else
    wlist=ws(:);
end
W=adjlist2mat(15,alist,wlist);
if isnorm
    L=l2nl(w2l(W));
else
    L=w2l(W);
end

[U0,D0]=eiga(L);
U=U0;
% get the "sparse" version of the GFT matrix
for i=1:14
    if abs(D0(i,i)-D0(i+1,i+1))<1e-4
        idx_1st_nz=find(all(abs(U(:,[i,i+1]))>1e-4,2),1);
        % apply a sparsifying rotation
        angle1=atan(U(idx_1st_nz,i+1)/U(idx_1st_nz,i));
        R2x2=[cos(angle1),-sin(angle1);sin(angle1),cos(angle1)];
        U(:,[i,i+1])=U(:,[i,i+1])*R2x2;
    end
end

% flip some signs 
for i=1:15
    idx_1st_nz=find(abs(U(:,i))>1e-4,1);
    if U(idx_1st_nz,i)<0
        U(:,i)=-U(:,i);
    end
end