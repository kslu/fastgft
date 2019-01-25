function [L,W]=unifgrid_general(N,wv,wh,wd,wa)
% create graph Laplacian matrix of an 8-connected grid with specified
% weights in vertical, horizontal, diagonal, and anti-diagonal directions, 
% respectively
% 
% L=unifgrid_general(N,wv,wh,wd,wa)
% 
% N: grid size
% wv: weights of vertical edges
% wh: weights of horizontal edges 
% wd: weights of diagonal edges
% wa: weights of anti-diagonal edges
% (wv, wh, wd, wa are all default to 1)
% 
% by: KS Lu
% 20180914
%

if nargin<5 || isempty(wa)
    wa=0;
end
if nargin<4 || isempty(wd)
    wd=0;
end
if nargin<3 || isempty(wh)
    wh=1;
end
if nargin<2 || isempty(wv)
    wv=1;
end

W=zeros(N^2);
W_lg=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
% horizontal and vertical edges
for i=1:N
    idx_ci=(i-1)*N+1:i*N;  % indices for column i
    idx_ri=i:N:N^2;  % indices for row i
    W(idx_ci,idx_ci)=W(idx_ci,idx_ci)+wv*W_lg;
    W(idx_ri,idx_ri)=W(idx_ri,idx_ri)+wh*W_lg;
end

% diagonal and anti-diagonal edges
for i=1:N-1
    % edges between the i-th and (i+1)-th columns, j-th and (j+1)-th rows
    for j=1:N-1
        idx_di=[(i-1)*N+j, i*N+j+1];
        idx_ai=[(i-1)*N+j+1, i*N+j];
        W(idx_di,idx_di)=W(idx_di,idx_di)+wd*[0,1;1,0];
        W(idx_ai,idx_ai)=W(idx_ai,idx_ai)+wa*[0,1;1,0];
    end
end

L=w2l(W);