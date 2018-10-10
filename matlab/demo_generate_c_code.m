%system('rm *.h');
addpath('third_party');
system('rm c_code/*');

%% skeletal graph 15
fprintf(sprintf('Generating code for sk15...\n'));

Wsk=skeleton15;
Lsk=w2l(Wsk);
[Usk,Dst]=eiga(Lsk);
Bsk=bmat(15,[4,5,6,10,11,12,15,14,13,9,8,7]);
Lsk_1=Bsk*Lsk*Bsk;
Lsk_p=Lsk_1([1,2,3,10,11,12,13,14,15],[1,2,3,10,11,12,13,14,15]);
Lsk_m=Lsk_1([4,5,6],[4,5,6]);  % Lm1=Lm2
[Usk_p,Dsk_p]=eiga(Lsk_p);
[Usk_m,Dsk_m]=eiga(Lsk_m);

print_c_gftmatrix(Usk,'sk15','c_code/sk15.txt');
print_c_gftmatrix(Usk_p/sqrt(2),'sk15_p','c_code/sk15_p.txt');
print_c_gftmatrix(Usk_m/sqrt(2),'sk15_m','c_code/sk15_m.txt');

%% skeletal graph 25
fprintf(sprintf('Generating code for sk25...\n'));

alist=[1,2,3,3,5,9,5,6,7,9,10,11,1,13,14,15,1,17,18,19,8,23,12,25;
       2,21,21,4,21,21,6,7,8,10,11,12,13,14,15,16,17,18,19,20,23,22,25,24]';
wlist=ones(1,24);
W=adjlist2mat(25,alist,wlist);
L=w2l(W);
[Usk25,Dsk25]=eiga(L);

idx_r=[5,6,7,8,22,23,13,14,15,16];
idx_l=[20,19,18,17,25,24,12,11,10,9];
idx_m=[1,2,3,4,21];
B1=bmat(25,[idx_l,idx_r]);
L1=B1*L*B1;

%figure; imagesc(L1([idx_m,idx_l,idx_r],[idx_m,idx_l,idx_r]));
idx_p=sort([idx_m, idx_l]);
idx_m1=[5,6,7,8,22,23];
idx_m2=[13,14,15,16];
Lp=L1(idx_p,idx_p);
Lm1=L1(idx_m1,idx_m1);
Lm2=L1(idx_m2,idx_m2);
[Usk25_p,Dsk25_p]=eiga(Lp);
[Usk25_m1,Dsk25_m1]=eiga(Lm1);
[Usk25_m2,Dsk25_m2]=eiga(Lm2);

print_c_gftmatrix(Usk25,'sk25','c_code/sk25.txt');
print_c_gftmatrix(Usk25_p/sqrt(2),'sk25_p','c_code/sk25_p.txt');
print_c_gftmatrix(Usk25_m1/sqrt(2),'sk25_m1','c_code/sk25_m1.txt');
print_c_gftmatrix(Usk25_m2/sqrt(2),'sk25_m2','c_code/sk25_m2.txt');

%% star graph 10
fprintf(sprintf('Generating code for star10...\n'));

n_star=10;
Wt=stargraph(n_star);
Lt=w2l(Wt);
[Ut,Dt]=eiga(Lt);
print_c_gftmatrix(Ut,'star10','c_code/star10.txt');

Wt_ppp=adjlist2mat(3,[1,1,1,3; 2,3,1,3]',[1,2*sqrt(2),8-2*sqrt(2),...
    1-2*sqrt(2)]);
Lt_ppp=w2l(Wt_ppp);
[Ut_ppp,Dt_ppp]=eiga(Lt_ppp);
print_c_gftmatrix(Ut_ppp/2/sqrt(2),'star10_ppp','c_code/star10_ppp.txt');

%% star graph 100
fprintf(sprintf('Generating code for star100...\n'));

W=stargraph(100);
L=w2l(W);
[Ustar100,Dstar100]=eiga(L);

B1=bmat(100,3:100);
B2=bmat(100,4:51);
B3=bmat(100,4:27);
B4=bmat(100,4:15);
B5=bmat(100,4:9);
B6=bmat(100,5:6);

L6=B6*B5*B4*B3*B2*B1*L*B1*B2*B3*B4*B5*B6;
%figure; imagesc(L6);

[Ustar100_pppppp,Dstar100_pppppp]=eiga(L6(1:5,1:5));
print_c_gftmatrix(Ustar100,'star100','c_code/star100.txt');
print_c_gftmatrix(Ustar100_pppppp/8,'star100_pppppp','c_code/star100_pppppp.txt');

%% bi-diagonally symmetric grid
fprintf(sprintf('Generating code for bd4x4...\n'));

a=0.5;
Lbd4x4=unifgrid_general(8,1,1,a,0);
[Ub,Db]=eiga(Lbd4x4);

Wbd4x4_pp=adjlist2mat(6,[1,2,3,1,2,2,3,5,1,2,3,4,5,6;2,3,4,5,5,6,6,6,...
    1,2,3,4,5,6]',[sqrt(2),1,sqrt(2),a,sqrt(2),sqrt(2)*a,sqrt(2),2,...
    2-sqrt(2),-(a+2)*(sqrt(2)-1),-2*(sqrt(2)-1),2-sqrt(2),2-sqrt(2),...
    (a+1)*(2-sqrt(2))]);
Wbd4x4_pm=adjlist2mat(4,[1,2,2,3,1,2,3,4;3,3,4,4,1,2,3,4]',[1,sqrt(2),...
    a,sqrt(2),2+2*a,4-sqrt(2)+2*a,2-2*sqrt(2)+a,2-sqrt(2)]);
Wbd4x4_mp=adjlist2mat(4,[1,1,2,2,1,2,3,4;2,3,3,4,1,2,3,4]',[1,sqrt(2)*a,...
    sqrt(2),sqrt(2),2-(sqrt(2)-1)*a,-2*(sqrt(2)-1),2+(2-sqrt(2))*(a+1),...
    2-sqrt(2)]);
Wbd4x4_mm=adjlist2mat(2,[1,1,2;2,1,2]',[1,2+2*a,2+a]);

Lbd4x4_pp=w2l(Wbd4x4_pp);
Lbd4x4_pm=w2l(Wbd4x4_pm);
Lbd4x4_mp=w2l(Wbd4x4_mp);
Lbd4x4_mm=w2l(Wbd4x4_mm);
[Ubd4x4_pp,Dbd4x4_pp]=eiga(Lbd4x4_pp);
[Ubd4x4_pm,Dbd4x4_pm]=eiga(Lbd4x4_pm);
[Ubd4x4_mp,Dbd4x4_mp]=eiga(Lbd4x4_mp);
[Ubd4x4_mm,Dbd4x4_mm]=eiga(Lbd4x4_mm);

print_c_gftmatrix(Ub,'bd4x4','c_code/bd4x4.txt');
print_c_gftmatrix(Ubd4x4_pp/2,'bd4x4_pp','c_code/bd4x4_pp.txt');
print_c_gftmatrix(Ubd4x4_pm/2,'bd4x4_pm','c_code/bd4x4_pm.txt');
print_c_gftmatrix(Ubd4x4_mp/2,'bd4x4_mp','c_code/bd4x4_mp.txt');
print_c_gftmatrix(Ubd4x4_mm/2,'bd4x4_mm','c_code/bd4x4_mm.txt');

%% 8x8 bidiagonally symmetric grid
fprintf(sprintf('Generating code for bd8x8...\n'));

L=unifgrid_general(8,1,1,a,0);
[Ubd8x8,Dbd8x8]=eiga(L);

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

L1=B2*B1*L*B1*B2;
%figure; imagesc(L1([idpp,idpm,idmp,idmm],[idpp,idpm,idmp,idmm]));

[Upp,Dpp]=eiga(L1(idpp,idpp));
[Upm,Dpm]=eiga(L1(idpm,idpm));
[Ump,Dmp]=eiga(L1(idmp,idmp));
[Umm,Dmm]=eiga(L1(idmm,idmm));

print_c_gftmatrix(Ubd8x8,'bd8x8','c_code/bd8x8.txt');
print_c_gftmatrix(Upp/2,'bd8x8_pp','c_code/bd8x8_pp.txt');
print_c_gftmatrix(Upm/2,'bd8x8_pm','c_code/bd8x8_pm.txt');
print_c_gftmatrix(Ump/2,'bd8x8_mp','c_code/bd8x8_mp.txt');
print_c_gftmatrix(Umm/2,'bd8x8_mm','c_code/bd8x8_mm.txt');

print_c_ptj(L,'bd8x8_ptj',[],[],[],120);
print_c_ptj(L1(idpp,idpp),'bd8x8_pp_ptj',[],[],[],100);
print_c_ptj(L1(idpm,idpm),'bd8x8_pm_ptj',[],[],[],100);
print_c_ptj(L1(idmp,idmp),'bd8x8_mp_ptj',[],[],[],100);
print_c_ptj(L1(idmm,idmm),'bd8x8_mm_ptj',[],[],[],100);
print_c_tj(L,'bd8x8_tj',[],[],[],4000);
print_c_tj(L1(idpp,idpp),'bd8x8_pp_tj',[],[],[],3000);
print_c_tj(L1(idpm,idpm),'bd8x8_pm_tj',[],[],[],3000);
print_c_tj(L1(idmp,idmp),'bd8x8_mp_tj',[],[],[],3000);
print_c_tj(L1(idmm,idmm),'bd8x8_mm_tj',[],[],[],3000);

%% 4x4 uniform grid
fprintf(sprintf('Generating code for dct4x4...\n'));

d1d=2-2*cos((0:3)*pi/4);  % eigenvalues for 1d
[xx,yy]=meshgrid(d1d,d1d);
d2d=xx(:)'+yy(:)';
Udct1d=dctmtx(4)';
Udct4x4=kron(Udct1d,Udct1d);
[~,d_idx]=sort(d2d);
Udct4x4=Udct4x4(:,d_idx);

Wdct4x4_pmp=adjlist2mat(2,[1,1;2,1]',[1,2]);
Wdct4x4_pmm=adjlist2mat(2,[1,1,2;2,1,2]',[1,4,2]);
Wdct4x4_mmp=adjlist2mat(3,[1,2,1,2,3;2,3,1,2,3]',[sqrt(2),sqrt(2),...
    6-sqrt(2),4-2*sqrt(2),2-sqrt(2)]);
Ldct4x4_pmp=w2l(Wdct4x4_pmp);
Ldct4x4_pmm=w2l(Wdct4x4_pmm);
Ldct4x4_mmp=w2l(Wdct4x4_mmp);
[Udct4x4_pmp,Ddct4x4_pmp]=eiga(Ldct4x4_pmp);
[Udct4x4_pmm,Ddct4x4_pmm]=eiga(Ldct4x4_pmm);
[Udct4x4_mmp,Ddct4x4_mmp]=eiga(Ldct4x4_mmp);
% Udct4x4_mpp is the same as Udct4x4_pmp
% Udct4x4_mpm is the same as Udct4x4_pmm

print_c_gftmatrix(Udct4x4,'dct4x4','c_code/dct4x4.txt');
print_c_gftmatrix(Udct4x4_pmp/2/sqrt(2),'dct4x4_pmp','c_code/dct4x4_pmp.txt');
print_c_gftmatrix(Udct4x4_pmm/2/sqrt(2),'dct4x4_pmm','c_code/dct4x4_pmm.txt');
print_c_gftmatrix(Udct4x4_mmp/2/sqrt(2),'dct4x4_mmp','c_code/dct4x4_mmp.txt');

%% 8x8 uniform grid
fprintf(sprintf('Generating code for dct8x8...\n'));

L=unifgrid_general(8);
%[Udct8x8,Ddct8x8]=eiga(L);

U1d=dctmtx(8)';
U2d=kron(U1d,U1d);
d1d=2-2*cos((0:7)*pi/8);  % eigenvalues for 1d
[xx,yy]=meshgrid(d1d,d1d);
d2d=xx(:)'+yy(:)';
[~,d_idx]=sort(d2d);
U2d=U2d(:,d_idx);
print_c_gftmatrix(U2d,'dct8x8','c_code/dct8x8.txt');

inv1=[1:32,40:-1:33,48:-1:41,56:-1:49,64:-1:57];
inv2=[1:4,9:12,17:20,25:28,33:36,41:44,49:52,57:60,...
    61:64,53:56,45:48,37:40,29:32,21:24,13:16,5:8];
inv3=[1:16,33,41,49,57,34,42,50,58,38,39,40,47,48,56,...
    63,62,54,61,53,45,59,51,43,35,60,52,44,36,24:-1:17,32:-1:25];
inv4=[1,2,9,10,17,18,25,26,5,6,7,8,33,41,49,57,...
    58,50,42,34,16:-1:13,27,28,19,20,11,12,3,4];
inv5=[1,2,3,4,17,25,20,27,26,18,12,11,10,9];
inv6=[1,9,10,2];

B1=bmat(64,inv1);
B2=bmat(64,inv2);
B3=bmat(64,inv3);
B4=bmat(64,inv4);
B5=bmat(64,inv5);
B6=bmat(64,inv6);

L6=B6*B5*B4*B3*B2*B1*L*B1*B2*B3*B4*B5*B6;

idx_new=[1,2,9,10,3,4,11,12,17,25,18,26,19,20,28,27,5:8,13:16,33:8:57,...
    34:8:58,21:24,29:32,35:8:59,36:8:60,37,38,39,40,46,47,48,55,56,64,...
    45,53,54,61,62,63];
%figure; imagesc(L6(idx_new,idx_new));

idx_pppmp=[3,4];
idx_pppmm=[11,12];
idx_ppmmp=[19,20,28];
idx_pmpp=[5,6,7,8];
idx_pmpm=[13,14,15,16];
idx_pmm=[21,22,23,24,29,30,31,32];
idx_mmp=[37,38,39,40,46,47,48,55,56,64];
idx_mmm=[45,53,54,61,62,63];

Ldct8x8_pppmp=L6(idx_pppmp,idx_pppmp);
Ldct8x8_pppmm=L6(idx_pppmm,idx_pppmm);
Ldct8x8_ppmpp=Ldct8x8_pppmp;
Ldct8x8_ppmpm=Ldct8x8_pppmm;
Ldct8x8_ppmmp=L6(idx_ppmmp,idx_ppmmp);
Ldct8x8_pmpp=L6(idx_pmpp,idx_pmpp);
Ldct8x8_pmpm=L6(idx_pmpm,idx_pmpm);
Ldct8x8_pmm=L6(idx_pmm,idx_pmm);
Ldct8x8_mppp=Ldct8x8_pmpp;
Ldct8x8_mppm=Ldct8x8_pmpm;
Ldct8x8_mpm=Ldct8x8_pmm;
Ldct8x8_mmp=L6(idx_mmp,idx_mmp);
Ldct8x8_mmm=L6(idx_mmm,idx_mmm);

[Udct8x8_pppmp,Ddct8x8_pppmp]=eiga(Ldct8x8_pppmp);
%[Udct8x8_pppmm,Ddct8x8_pppmm]=eiga(Ldct8x8_pppmm); % the same
[Udct8x8_ppmmp,Ddct8x8_ppmmp]=eiga(Ldct8x8_ppmmp);
[Udct8x8_pmpp,Ddct8x8_pmpp]=eiga(Ldct8x8_pmpp);
[Udct8x8_pmpm,Ddct8x8_pmpm]=eiga(Ldct8x8_pmpm);
[Udct8x8_pmm,Ddct8x8_pmm]=eiga(Ldct8x8_pmm);
[Udct8x8_mmp,Ddct8x8_mmp]=eiga(Ldct8x8_mmp);
[Udct8x8_mmm,Ddct8x8_mmm]=eiga(Ldct8x8_mmm);

print_c_gftmatrix(Udct8x8_pppmp/4/sqrt(2),'dct8x8_pppmp','c_code/dct8x8_pppmp.txt');
%print_c_gftmatrix(Udct8x8_pppmm/4/sqrt(2),'dct8x8_pppmm','c_code/dct8x8_pppmm.txt');
print_c_gftmatrix(Udct8x8_ppmmp/4/sqrt(2),'dct8x8_ppmmp','c_code/dct8x8_ppmmp.txt');
print_c_gftmatrix(Udct8x8_pmpp/4,'dct8x8_pmpp','c_code/dct8x8_pmpp.txt');
print_c_gftmatrix(Udct8x8_pmpm/4,'dct8x8_pmpm','c_code/dct8x8_pmpm.txt');
print_c_gftmatrix(Udct8x8_pmm/2/sqrt(2),'dct8x8_pmm','c_code/dct8x8_pmm.txt');
print_c_gftmatrix(Udct8x8_mmp/2/sqrt(2),'dct8x8_mmp','c_code/dct8x8_mmp.txt');
print_c_gftmatrix(Udct8x8_mmm/2/sqrt(2),'dct8x8_mmm','c_code/dct8x8_mmm.txt');

%% 4x4 Z-shaped graph
fprintf(sprintf('Generating code for z4x4...\n'));

wd=2;
Lz4=unifgrid_general(4,0,1,0,wd);
[Uz4,Dz4]=eiga(Lz4);
B16=bmat(16);
Lz4_1=B16*Lz4*B16;
[Uz4_p,Dz4_p]=eiga(Lz4_1(1:8,1:8));
[Uz4_m,Dz4_m]=eiga(Lz4_1(9:16,9:16));

print_c_gftmatrix(Uz4,'z4x4','c_code/z4x4.txt');
print_c_gftmatrix(Uz4_p,'z4x4_p','c_code/z4x4_p.txt');
print_c_gftmatrix(Uz4_m,'z4x4_m','c_code/z4x4_m.txt');

%% 8x8 Z-shaped graph
fprintf(sprintf('Generating code for z8x8...\n'));
Lz8=unifgrid_general(8,0,1,0,wd);
[Uz8,Dz8]=eiga(Lz8);
B64=bmat(64);
Lz8_1=B64*Lz8*B64;
[Uz8_p,Dz8_p]=eiga(Lz8_1(1:32,1:32));
[Uz8_m,Dz8_m]=eiga(Lz8_1(33:64,33:64));

print_c_gftmatrix(Uz8,'z8x8','c_code/z8x8.txt');
print_c_gftmatrix(Uz8_p,'z8x8_p','c_code/z8x8_p.txt');
print_c_gftmatrix(Uz8_m,'z8x8_m','c_code/z8x8_m.txt');

[~,z4x4_isortd]=sort([diag(Dz4_p)',diag(Dz4_m)']);
oidx_z4x4=zeros(1,16); oidx_z4x4(z4x4_isortd)=1:16;
[~,z8x8_isortd]=sort([diag(Dz8_p)',diag(Dz8_m)']);
oidx_z8x8=zeros(1,64); oidx_z8x8(z8x8_isortd)=1:64;

print_c_ptj(Lz8,'z8x8_ptj',[],[],[],120);
print_c_ptj(Lz8_1(1:32,1:32),'z8x8_p_ptj',[],[],[],100);
print_c_ptj(Lz8_1(33:64,33:64),'z8x8_m_ptj',[],[],[],100);
print_c_tj(Lz8,'z8x8_tj',[],[],[],4000);
print_c_tj(Lz8_1(1:32,1:32),'z8x8_p_tj',[],[],[],3000);
print_c_tj(Lz8_1(33:64,33:64),'z8x8_m_tj',[],[],[],3000);


%% move the c_code folder
system('mv *tj_coords.h c_code/');
system('mv *tj_angles.h c_code/');
system('mv *tj_idx.h c_code/');

%% concatenate all c code
fprintf(sprintf('Combining files...\n'));

files_mat_cat={'star10','star10_ppp','star100','star100_pppppp',...
    'bd4x4','bd4x4_pp','bd4x4_pm','bd4x4_mp','bd4x4_mm',...
    'bd8x8','bd8x8_pp','bd8x8_pm','bd8x8_mp','bd8x8_mm',...
    'dct4x4','dct4x4_pmp','dct4x4_mmp','dct8x8','dct8x8_pppmp',...
    'dct8x8_ppmmp','dct8x8_pmpp','dct8x8_pmm','dct8x8_mmp','dct8x8_mmm',...
    'sk15','sk15_p','sk15_m','sk25','sk25_p','sk25_m1','sk25_m2',...
    'z4x4','z4x4_p','z4x4_m','z8x8','z8x8_p','z8x8_m'};

out_mat='c_code/mat_all.h';
for i=1:length(files_mat_cat)
    fstr_mat=sprintf('c_code/%s.txt',files_mat_cat{i});
    system(['cat ', fstr_mat, ' >> ', out_mat]);
end

files_tj_cat={'bd8x8','bd8x8_pp','bd8x8_pm','bd8x8_mp','bd8x8_mm','z8x8',...
    'z8x8_p','z8x8_m'};

out_tj='c_code/tj_all.h';
for i=1:length(files_tj_cat)
    fstr_tj_c=sprintf('c_code/%s_tj_coords.h',files_tj_cat{i});
    fstr_tj_a=sprintf('c_code/%s_tj_angles.h',files_tj_cat{i});
    fstr_tj_i=sprintf('c_code/%s_tj_idx.h',files_tj_cat{i});
    system(['cat ', fstr_tj_c, ' >> ', out_tj]);
    system(['cat ', fstr_tj_a, ' >> ', out_tj]);
    system(['cat ', fstr_tj_i, ' >> ', out_tj]);
end

out_ptj='c_code/ptj_all.h';
for i=1:length(files_tj_cat)
    fstr_ptj_c=sprintf('c_code/%s_ptj_coords.h',files_tj_cat{i});
    fstr_ptj_a=sprintf('c_code/%s_ptj_angles.h',files_tj_cat{i});
    fstr_ptj_i=sprintf('c_code/%s_ptj_idx.h',files_tj_cat{i});
    system(['cat ', fstr_ptj_c, ' >> ', out_ptj]);
    system(['cat ', fstr_ptj_a, ' >> ', out_ptj]);
    system(['cat ', fstr_ptj_i, ' >> ', out_ptj]);
end

% Then, run this in command line:
% cp ~/Dropbox/MATLAB/STAC/symdecomp/c_code/tj_all.h givens_tj.h && cp ~/Dropbox/MATLAB/STAC/symdecomp/c_code/ptj_all.h givens_ptj.h