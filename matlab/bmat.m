function B=bmat(N,ordered_idx)
% B matrix (a scaled butterfly)
% 
% B=bmat(N,ordered_idx)
% 
% N: size, must be even
% B: butterfly matrix
% ordered_idx (optional): indicate the pairing scenario. It can include
%     permutations and node selection. For example when N=6,
%     ordered_idx=[1,3,4,6], it means that node 1 and 6 are paired, node 3
%     and 4 are paired, and processed using a butterfly, and no operation 
%     is applied to nodes 2 and 5. 
% 

if nargin>1
    if size(ordered_idx,1)~=1 && size(ordered_idx,2)>2
        error('ordered_idx cannot have more than two columns');
    elseif size(ordered_idx,2)==2
        ordered_idx=[ordered_idx(:,1);ordered_idx(end:-1:1,2)];
    end
    Bcomp=bmat(numel(ordered_idx));
    B=eye(N);
    B(ordered_idx,ordered_idx)=Bcomp;
else
    if mod(N,2)==1
        h=(N-1)/2;
        B=1/sqrt(2)*[    eye(h), zeros(h,1),    jmat(h); 
                     zeros(1,h),    sqrt(2), zeros(1,h);
                        jmat(h), zeros(h,1),    -eye(h)];
    else
        B=1/sqrt(2)*[eye(N/2), jmat(N/2); jmat(N/2), -eye(N/2)];
    end
end