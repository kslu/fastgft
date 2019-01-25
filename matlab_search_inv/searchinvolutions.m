function valid_inv=searchinvolutions(W)
% SEARCHINVOLUTIONS searches valid involutions of a graph using degree lists
% 
% valid_inv=searchinvolutions(W)
% 
% 20190124
%
addpath('third_party/');
n=size(W,1);
maxinvs=500000;

% degree lists
dv=sum(W);
fv=sum(W~=0);

%%
[C,ia,ic]=unique([dv;fv]','rows');
nbin=size(C,1);  % number of bins of the histogram
bins=cell(1,nbin);
binsizes=zeros(1,nbin);
for i=1:nbin
    binsizes(i)=nnz(ic==i);
    bins{i}=find(ic==i);
end

Ts=[1, 2, 4, 10, 26, 76, 232, 764, 2620, 9496, 35696, 140152];
ninvs2search=prod(Ts(binsizes));
if ninvs2search>maxinvs
    error('Search dimension too large!');
end

nbin2=nnz(binsizes>1);
N4search=bins(binsizes>1);
I4search=cell(1,nbin2);
idx4search=cell(1,nbin2);
for i=1:length(N4search)
    temp=allinvolutions(numel(N4search{i}));
    idx_bin_i=N4search{i};
    I4search{i}=idx_bin_i(temp);
    idx4search{i}=1:size(temp,1);
end

idx_all=allcomb2(idx4search);

%% search
inv0=1:n;
valid_inv=[];
for i=1:ninvs2search
    vec_idx=idx_all(i,:);
    
    for j=1:nbin2
        invs_j=I4search{j};
        inv0(invs_j(1,:))=invs_j(vec_idx(j),:);
    end
    
    if all(all(W==W(inv0,inv0)))
        valid_inv=[valid_inv; inv0];
    end
end
