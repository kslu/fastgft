function [W,L]=stargraph(n)
W=zeros(n);
W(1,2:end)=1;
W(2:end,1)=1;
L=w2l(W);