function [alist,wlist,n]=mat2adjlist(W)

n=size(W,1);
wlist=W(find(triu(W)));
[ii,jj]=find(triu(W));
alist=[ii,jj];