function allinvs=allinvolutions(n)

if n==1
    allinvs=1;
elseif n==2
    allinvs=[1,2;2,1];
else
    % recursive formula: T(n)=T(n-1)+(n-1)T(n-2)
    %   The first term: all involutions with the last element fixed
    %   The second term: all involutions with the last element moved

    % recursive implementation
%     A1=allinvolutions(n-1);
%     A2=allinvolutions(n-2);
%     T1=size(A1,1);
%     T2=size(A2,1);
% 
%     allinvs=zeros(T1+(n-1)*T2,n);
%     allinvs(1:T1,1:n-1)=A1;
%     allinvs(1:T1,end)=n;
%     for i=1:n-1
%         idxr=T1+(i-1)*T2+1:T1+i*T2;
%         allinvs(idxr,i)=n;
%         allinvs(idxr,n)=i;
%         idxc_res=setdiff(1:n,[i,n]);
%         idx_no_i=setdiff(1:n,i);
%         allinvs(idxr,idxc_res)=idx_no_i(A2);
%     end
    
    % dynamic programming implementation
    Ajm1=[1,2;2,1];
    Ajm2=1;
    Tjm1=2;
    Tjm2=1;
    for j=3:n
        [allinvs,T]=recursion_helper(Ajm1,Ajm2,Tjm1,Tjm2,j);
        Ajm2=Ajm1;
        Ajm1=allinvs;
        Tjm2=Tjm1;
        Tjm1=T;
    end
end

function [A,T]=recursion_helper(A1,A2,T1,T2,n)
T=T1+(n-1)*T2;
A=zeros(T,n);
A(1:T1,1:n-1)=A1;
A(1:T1,end)=n;
for i=1:n-1
    idxr=T1+(i-1)*T2+1:T1+i*T2;
    A(idxr,i)=n;
    A(idxr,n)=i;
    idxc_res=setdiff(1:n,[i,n]);
    idx_no_i=setdiff(1:n,i);
    A(idxr,idxc_res)=idx_no_i(A2);
end