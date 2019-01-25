% example of the tree involution search in the report
n=13;
alist=[1,2;2,3;3,4;3,5;2,6;1,7;7,8;7,9;9,10;9,11;1,12;12,13];
wlist=[2,3,1,1,1,2,1,3,1,1,1,1];
W=adjlist2mat(13,alist,wlist);

% using tree search
[invs,ninv,desc]=treeinvolution(13,W);

% using degree lists
invs_d=searchinvolutions(W);