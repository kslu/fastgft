% for: get the relative error and number of multiplications for each GFT 
% implementation
close all;
clear;

n=64;
num_layers=120;
num_layers_subgft=100;

re_b_ptj_per_layer=zeros(1,num_layers+1);
re_b_btfptj_per_layer=zeros(1,num_layers_subgft+1);
re_z_ptj_per_layer=zeros(1,num_layers+1);
re_z_btfptj_per_layer=zeros(1,num_layers_subgft+1);
mult_per_layer=zeros(1,num_layers);
mult_per_layer_subgft=zeros(1,num_layers_subgft);

%% A. bi-diagonally symmetric 8x8 grid
a=0.5;
Lb=unifgrid_general(8,1,1,a,0);
[Ub,Db]=eiga(Lb);

idll=[2:8,11:16,20:24,29:32,38:40,47,48,56];
idur=[9:8:57,18:8:58,27:8:59,36:8:60,45,53,61,54,62,63];
idul=[1:7,9:14,17:21,25:28,33:35,41,42,49];
idlr=[64:-8:16,63:-8:23,62:-8:30,61:-8:37,60,52,44,59,51,58];
B1=bmat(64,[idll, fliplr(idur)]);
B2=bmat(64,[idul, fliplr(idlr)]);
idpp=[1:8,10:15,19:22,28,29];
idpm=[16,23,24,30,31,32,37,38,39,40,46,47,48,55,56,64];
idmp=[9,17,18,25,26,27,33,34,35,36,41,42,43,49,50,57];
idmm=[44,45,51:54,58:63];
L1=B2*B1*Lb*B1*B2;
Lbpp=L1(idpp,idpp);
Lbpm=L1(idpm,idpm);
Lbmp=L1(idmp,idmp);
Lbmm=L1(idmm,idmm);
[Upp,Dpp]=eiga(Lbpp);
[Upm,Dpm]=eiga(Lbpm);
[Ump,Dmp]=eiga(Lbmp);
[Umm,Dmm]=eiga(Lbmm);

% A.1 Run parallel truncated Jacobi algorithm
[facs_b,~,~,~,~,idx_b]=givens_ptj(Lb,num_layers*n/2,n/2);
[~,idx_b_orig]=sort(diag(Lb)');  idx_b=[idx_b_orig; idx_b];

[facs_b_pp,~,~,~,~,idx_b_pp]=...
    givens_ptj(Lbpp,num_layers_subgft*n/2,n/2);
[~,idx_b_pp_orig]=sort(diag(Lbpp)');  idx_b_pp=[idx_b_pp_orig; idx_b_pp];
[facs_b_pm,~,~,~,~,idx_b_pm]=...
    givens_ptj(Lbpm,num_layers_subgft*n/2,n/2);
[~,idx_b_pm_orig]=sort(diag(Lbpm)');  idx_b_pm=[idx_b_pm_orig; idx_b_pm];
[facs_b_mp,~,~,~,~,idx_b_mp]=...
    givens_ptj(Lbmp,num_layers_subgft*n/2,n/2);
[~,idx_b_mp_orig]=sort(diag(Lbmp)');  idx_b_mp=[idx_b_mp_orig; idx_b_mp];
[facs_b_mm,~,~,~,~,idx_b_mm]=...
    givens_ptj(Lbmm,num_layers_subgft*n/2,n/2);
[~,idx_b_mm_orig]=sort(diag(Lbmm)');  idx_b_mm=[idx_b_mm_orig; idx_b_mm];

% A.2 Find the numbers of multiplications
% A.2.1 overall GFT
for i=1:num_layers
    mult_per_layer(i)=4*nnz(triu(facs_b{i},1));
end
nmult_bd8x8_ptj=[0,cumsum(mult_per_layer)];
nmult_bd8x8_mat=nnz(Ub);
nmult_bd8x8_btf=nnz(Upp)+nnz(Upm)+nnz(Ump)+nnz(Umm);
% A.2.2 sub-GFTs
for i=1:num_layers_subgft
    mult_per_layer_subgft(i)=4*nnz(triu(facs_b_pp{i},1));
end
nmult_bd8x8_ptj_pp=[0,cumsum(mult_per_layer_subgft)];
for i=1:num_layers_subgft
    mult_per_layer_subgft(i)=4*nnz(triu(facs_b_pm{i},1));
end
nmult_bd8x8_ptj_pm=[0,cumsum(mult_per_layer_subgft)];
for i=1:num_layers_subgft
    mult_per_layer_subgft(i)=4*nnz(triu(facs_b_mp{i},1));
end
nmult_bd8x8_ptj_mp=[0,cumsum(mult_per_layer_subgft)];
for i=1:num_layers_subgft
    mult_per_layer_subgft(i)=4*nnz(triu(facs_b_mm{i},1));
end
nmult_bd8x8_ptj_mm=[0,cumsum(mult_per_layer_subgft)];

% A.3 Find the REs (get equivalent GFT matrices first)
Ub_ptj=eye(n);
Uzp_btfptj=eye(20);
Uzm_btfptj=eye(16);
Ubmp_btfptj=eye(16);
Ubmm_btfptj=eye(12);
% A.3.1 full GFT
re_b_ptj_per_layer(1)=norm(abs(Ub)-eye(n),'fro')/sqrt(n);  % zero layers
for i=1:num_layers
    Ub_ptj=Ub_ptj*full(facs_b{i});
    re_b_ptj_per_layer(i+1)=norm(abs(Ub_ptj(:,idx_b(i+1,:))'*Ub)-eye(n),'fro')/sqrt(n);
    %re_b_ptj_per_layer(i+1)=norm(abs(Ub_ptj'*Ub)-eye(n),'fro')/sqrt(n);
end
% A.3.2 sub-GFT
re_b_btfptj_per_layer(1)=sqrt(norm(abs(Upp)-eye(20),'fro')^2+...
    norm(abs(Upm)-eye(16),'fro')^2+norm(abs(Ump)-eye(16),'fro')^2+...
    norm(abs(Umm)-eye(12),'fro')^2)/sqrt(n);  % zero layers
for i=1:num_layers_subgft
    Uzp_btfptj=Uzp_btfptj*full(facs_b_pp{i});
    Uzm_btfptj=Uzm_btfptj*full(facs_b_pm{i});
    Ubmp_btfptj=Ubmp_btfptj*full(facs_b_mp{i});
    Ubmm_btfptj=Ubmm_btfptj*full(facs_b_mm{i});
    re_b_btfptj_per_layer(i+1)=sqrt(...
        norm(abs(Uzp_btfptj(:,idx_b_pp(i+1,:))'*Upp)-eye(20),'fro')^2+...
        norm(abs(Uzm_btfptj(:,idx_b_pm(i+1,:))'*Upm)-eye(16),'fro')^2+...
        norm(abs(Ubmp_btfptj(:,idx_b_mp(i+1,:))'*Ump)-eye(16),'fro')^2+...
        norm(abs(Ubmm_btfptj(:,idx_b_mm(i+1,:))'*Umm)-eye(12),'fro')^2)/sqrt(n);
end


%% B. Z-shaped 8x8 grid
wd=2;
Lz=unifgrid_general(8,0,1,0,wd);
[Uz,Dz]=eiga(Lz);
B64=bmat(64);
Lz_1=B64*Lz*B64;
Lzp=Lz_1(1:32,1:32);
Lzm=Lz_1(33:64,33:64);
[Uzp,Dzp]=eiga(Lzp);
[Uzm,Dzm]=eiga(Lzm);

% B.1 Run parallel truncated Jacobi algorithm
[facs_z,~,~,~,~,idx_z]=givens_ptj(Lz,num_layers*n/2,n/2);
[~,idx_z_orig]=sort(diag(Lz)');  idx_z=[idx_z_orig; idx_z];
[facs_z_p,~,~,~,~,idx_z_p]=givens_ptj(Lzp,num_layers_subgft*n/2,n/2);
[~,idx_z_p_orig]=sort(diag(Lzp)');  idx_z_p=[idx_z_p_orig; idx_z_p];
[facs_z_m,~,~,~,~,idx_z_m]=givens_ptj(Lzm,num_layers_subgft*n/2,n/2);
[~,idx_z_m_orig]=sort(diag(Lzm)');  idx_z_m=[idx_z_m_orig; idx_z_m];

% B.2 Find the numbers of multiplications
% B.2.1 overall GFT
for i=1:num_layers
    mult_per_layer(i)=4*nnz(triu(facs_z{i},1));
end
nmult_z8x8_ptj=[0,cumsum(mult_per_layer)];
nmult_z8x8_mat=nnz(Uz);
nmult_z8x8_btf=nnz(Uzp)+nnz(Uzm);
% B.2.2 sub-GFTs
for i=1:num_layers_subgft
    mult_per_layer_subgft(i)=4*nnz(triu(facs_z_p{i},1));
end
nmult_z8x8_ptj_p=[0,cumsum(mult_per_layer_subgft)];
for i=1:num_layers_subgft
    mult_per_layer_subgft(i)=4*nnz(triu(facs_z_m{i},1));
end
nmult_z8x8_ptj_m=[0,cumsum(mult_per_layer_subgft)];

% B.3 Find the REs (get equivalent GFT matrices first)
Uz_ptj=eye(n);
Uzp_btfptj=eye(32);
Uzm_btfptj=eye(32);
% B.3.1 full GFT
re_z_ptj_per_layer(1)=norm(abs(Uz)-eye(n),'fro')/sqrt(n);  % zero layers
for i=1:num_layers
    Uz_ptj=Uz_ptj*full(facs_z{i});
    re_z_ptj_per_layer(i+1)=norm(abs(Uz_ptj(:,idx_z(i+1,:))'*Uz)-eye(n),'fro')/sqrt(n);
    %re_z_ptj_per_layer(i+1)=norm(abs(Uz_ptj'*Uz)-eye(n),'fro')/sqrt(n);
end
% B.3.2 sub-GFTs
re_z_btfptj_per_layer(1)=sqrt(norm(abs(Uzp)-eye(32),'fro')^2+...
    norm(abs(Uzm)-eye(32),'fro')^2)/sqrt(n);  % zero layers
for i=1:num_layers_subgft
    Uzp_btfptj=Uzp_btfptj*full(facs_z_p{i});
    Uzm_btfptj=Uzm_btfptj*full(facs_z_m{i});
    re_z_btfptj_per_layer(i+1)=sqrt(...
        norm(abs(Uzp_btfptj(:,idx_z_p(i+1,:))'*Uzp)-eye(32),'fro')^2+...
        norm(abs(Uzm_btfptj(:,idx_z_m(i+1,:))'*Uzm)-eye(32),'fro')^2)/sqrt(n);
end

%% save the results
save('nmults','nmult_bd8x8_mat','nmult_bd8x8_btf','nmult_bd8x8_ptj',...
    'nmult_bd8x8_ptj_pp','nmult_bd8x8_ptj_pm','nmult_bd8x8_ptj_mp',...
    'nmult_bd8x8_ptj_mm','nmult_z8x8_mat','nmult_z8x8_btf',...
    'nmult_z8x8_ptj','nmult_z8x8_ptj_p','nmult_z8x8_ptj_m');
save('res','re_b_ptj_per_layer','re_b_btfptj_per_layer',...
    're_z_ptj_per_layer','re_z_btfptj_per_layer');