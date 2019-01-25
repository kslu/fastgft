function [ilist,ninv,blist]=treeinvolution(n,W,verbose)
% TREEINVOLUTION returns all involutions phi that gives a graph symmetry
% on a tree
% 
% [ilist,ninv,blist]=treeinvolution(n,W,verbose)
% 
% n: number of nodes
% W: adjacency matrix of a tree graph
% ilist: list of all valid involutions
% ninv: number of involutions
% blist: list of branch descriptors
% 
% 20190124

if nargin<3
    verbose=0;
end
rt=1;
alist=mat2adjlist(W);
Gt=graph(alist(:,1),alist(:,2));

%% find the center(s)
Tb1=bfsearch(Gt,rt);
Tb2=bfsearch(Gt,Tb1(end),'edgetonew');
% list of parents
plist=zeros(1,n);
plist(Tb2(:,2))=Tb2(:,1);
% find longest path: start with the last node, find parents until reaching
% the first node
n_cur=Tb2(end,2);  % last node (of Tb2)
maxpath=[n_cur, zeros(1,n-1)];
len_diameter=1;
while plist(n_cur)~=0
    len_diameter=len_diameter+1;
    maxpath(len_diameter)=plist(n_cur);
    n_cur=plist(n_cur);
end
maxpath=maxpath(1:len_diameter);
if mod(len_diameter,2)==1
    centers=maxpath((len_diameter+1)/2);
else
    centers=maxpath((len_diameter+[0,2])/2);
end

%% bfs again on the rooted tree, and create descriptors bottom-up
Tc=bfsearch(Gt,centers(1),'edgetonew');
bfslist=[centers(1),Tc(:,2)'];
plist=zeros(1,n);  % list of parents
plist(Tc(:,2))=Tc(:,1);
clist=cell(1,n);  % list of children
blist=cell(1,n);  % list of branch descriptors
tlist=cell(1,n);  % list of nodes in the outgoing branch
%ilist=cell(1,n);  % list of valid involutions
ninv=0;
for i=n:-1:1
    n_cur=bfslist(i);
    cur_clist=clist{n_cur};
    if ~isempty(cur_clist)
        n_child=numel(cur_clist);
        templist=cell(1,n_child);
        for j=1:numel(cur_clist)
            wi=W(n_cur,cur_clist(j));
            templist{j}=['(',sprintf('%d',wi),blist{cur_clist(j)},')'];
        end
        [sorted_branches,idx_childs]=sort(templist);
        % check identical branches (sort strings, and then find 
        % identical consecutive strings)
        childlist=clist{n_cur};
        sorted_childlist=childlist(idx_childs);
        for k=1:n_child-1
            if strcmp(sorted_branches(k),sorted_branches(k+1))
                ninv=ninv+1;
                ilist{ninv}=[sorted_childlist(k), tlist{sorted_childlist(k)};
                             sorted_childlist(k+1), tlist{sorted_childlist(k+1)}];
            end
        end
        % concatenate childlists 
        cur_tlist=zeros(1,n);
        iel=1;
        for l=1:n_child
            cur_tl=tlist{sorted_childlist(l)};
            len_cur_tl=numel(cur_tl);
            cur_tlist(iel:iel+len_cur_tl)=[sorted_childlist(l), cur_tl];
            iel=iel+len_cur_tl+1;
        end
        cur_tlist=cur_tlist(1:iel-1);
        tlist{n_cur}=cur_tlist;
        % concatenate branch descriptors
        blist{n_cur}=strjoin(sorted_branches,'');
    else
        blist{n_cur}='';
    end
    % fill its parent's children list
    if i~=1
        clist{plist(n_cur)}=[clist{plist(n_cur)}, n_cur];
    end
    
    if verbose
        fprintf('.');
        if mod(n+1-i,50)==0, fprintf('\n'); end
    end
end

%% check the other center for complete symmetry
% For centers c1, c2, assume the node on the other side of c1 is b1, 
% we check if they are equal:
% branch 1: c1 -> b1 -> ...
% branch 2: c2 -> ...
if numel(centers)==2
    c1=centers(1);
    c2=centers(2);
    clist1=clist{c1};
    b1=clist1(clist1~=c2);
    wb=W(c1,b1);
    % extract two branch descriptors
    blist_c2=blist{c2};
    blist_b1=['(', sprintf('%d',wb), blist{b1}, ')'];
    if strcmp(blist_c2,blist_b1)
        ninv=ninv+1;
        ilist{ninv}=[c1, b1, tlist{b1}; c2, tlist{c2};];
    end
end

if ninv==0
    ilist=[];
end