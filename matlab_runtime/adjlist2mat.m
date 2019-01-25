function A = adjlist2mat(n,alist,wlist,is_directed)
% ADJLIST2MAT converts an adjacency list to a weight matrix
% 
% A = adjlist2mat(n,alist,wlist,is_directed)
% 
% n: number of vertices
% alist: an mx2 adjacency list matrix. m is the number of edges
% wlist: a mx1 vector of weights
% 
% 20180419
if nargin<4
    is_directed=0;
end
if size(alist,1)==2 && size(alist,2)~=2
    alist=alist';
end
m=size(alist,1);
if nargin<3
    wlist=ones(m,1);
end
A=zeros(n);
src2dst = sub2ind([n,n],alist(:,1),alist(:,2));
dst2src = sub2ind([n,n],alist(:,2),alist(:,1));
A(src2dst) = wlist(:);
if ~is_directed
    A(dst2src) = wlist(:);
end