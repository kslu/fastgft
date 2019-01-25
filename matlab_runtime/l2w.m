function W=l2w(L)
% L2W Laplacian matrix to weight matrix
% 
% W=l2w(L)
% 
% by: KS Lu
% 20170302
%
if ~islaplacian(L)
    %error('L is not a valid Laplacian matrix');
end
NN=size(L,1);
W=-L;
W(1:NN+1:end)=sum(L);