% for: plot runtime vs error of different fast GFT algorithms
close all;
clear;

dirpath='~/Dropbox/mycode/fastgft/output/';
files={'bd8x8', 'z8x8'};
load('nmults.mat');
load('res.mat');
%files={'bd8x8'};

xrange_mult=7600;
xrange_rt=3e-5;

for i=1:length(files)
    % first pass: get length of each fast GFT method
    fid=fopen([dirpath, files{i}, '.txt']);
    for j=1:3, fgetl(fid); end
    nrow_tj=0;
    nrow_ptj=0;
    nrow_btftj=0;
    nrow_btfptj=0;
    while 1
        li=fgetl(fid);
        if li==-1, break; end
        li_list=strsplit(li,{':'});
        switch li_list{1}
            case 'TJ-GFT'
                nrow_tj=nrow_tj+1;
            case 'PTJ-GFT'
                nrow_ptj=nrow_ptj+1;
            case 'BTFTJ-GFT'
                nrow_btftj=nrow_btftj+1;
            case 'BTFPTJ-GFT'
                nrow_btfptj=nrow_btfptj+1;
        end
    end
    
    % pass 2: collect results
    rt_tj=zeros(1,nrow_tj);
    err_tj=zeros(1,nrow_tj);
    ngivens_tj=cell(1,nrow_tj);
    rt_ptj=zeros(1,nrow_ptj);
    err_ptj=zeros(1,nrow_ptj);
    nlayers_ptj=cell(1,nrow_ptj);
    rt_btftj=zeros(1,nrow_btftj);
    err_btftj=zeros(1,nrow_btftj);
    ngivens_btftj=cell(1,nrow_btftj);
    rt_btfptj=zeros(1,nrow_btfptj);
    err_btfptj=zeros(1,nrow_btfptj);
    nlayers_btfptj=cell(1,nrow_btfptj);
    
    fid=fopen([dirpath, files{i}, '.txt']);
    % line 1: number of inputs
    l1_list=strsplit(fgetl(fid),{' ','='});
    n_input=str2double(l1_list{2});
    % line 2: matrix GFT
    l2_list=strsplit(fgetl(fid),{' ',':'});
    rt_mat=str2double(l2_list{3});
    % line 3: butterfly + matrix GFT
    l3_list=strsplit(fgetl(fid),{' ',':'});
    rt_btf=str2double(l3_list{3});
    % TJ-GFT
    for k=1:nrow_tj
        tj_list=strsplit(fgetl(fid),{' ',':',',','=','(',')'});
        rt_tj(k)=str2double(tj_list{2});
        err_tj(k)=str2double(tj_list{7});
        ngivens_tj{k}=tj_list{3};
    end
    % PTJ-GFT
    for k=1:nrow_ptj
        ptj_list=strsplit(fgetl(fid),{' ',':',',','=','(',')'});
        rt_ptj(k)=str2double(ptj_list{2});
        err_ptj(k)=str2double(ptj_list{6});
        nlayers_ptj{k}=ptj_list{3};
    end
    % BTF-TJ-GFT
    for k=1:nrow_btftj
        btftj_list=strsplit(fgetl(fid),{' ',':',',','=','(',')'});
        rt_btftj(k)=str2double(btftj_list{2});
        err_btftj(k)=str2double(btftj_list{7});
        ngivens_btftj{k}=btftj_list{3};
    end
    % BTF-PTJ-GFT
    for k=1:nrow_btfptj
        btfptj_list=strsplit(fgetl(fid),{' ',':',',','=','(',')'});
        rt_btfptj(k)=str2double(btfptj_list{2});
        err_btfptj(k)=str2double(btfptj_list{6});
        nlayers_btfptj{k}=btfptj_list{3};
    end
    
    % range control (do not show text out of this range)
    idx_show_ptj=rt_ptj/n_input<xrange_rt;
    idx_show_btfptj=rt_btfptj/n_input<xrange_rt;
    
    %% plot #multiplication vs error
    switch files{i}
        case 'bd8x8'
            nmult_mat=nmult_bd8x8_mat(1:5:end);
            nmult_btf=nmult_bd8x8_btf(1:5:end);
            nmult_ptj=nmult_bd8x8_ptj(1:5:end);
            nmult_btfptj=nmult_bd8x8_ptj_pp(1:5:end)+...
                nmult_bd8x8_ptj_pm(1:5:end)+nmult_bd8x8_ptj_mp(1:5:end)+...
                nmult_bd8x8_ptj_mm(1:5:end);
        case 'z8x8'
            nmult_mat=nmult_z8x8_mat(1:5:end);
            nmult_btf=nmult_z8x8_btf(1:5:end);
            nmult_ptj=nmult_z8x8_ptj(1:5:end);
            nmult_btfptj=nmult_z8x8_ptj_p(1:5:end)+nmult_z8x8_ptj_m(1:5:end);
    end
    % range control (do not show text out of this range)
    idx_show_mult_ptj=nmult_ptj<xrange_mult;
    idx_show_mult_btfptj=nmult_btfptj<xrange_mult;
    
    figure;
    plot(nmult_mat,0,'*','linewidth',4,'markersize',20); 
    hold on; grid on;
    plot(nmult_btf,0,'s','linewidth',4,'markersize',20);
    plot(nmult_ptj,err_ptj,'+:','linewidth',2,'markersize',10);
    text(nmult_ptj(idx_show_mult_ptj)+50,err_ptj(idx_show_mult_ptj)+...
        0.8,nlayers_ptj(idx_show_mult_ptj),'fontsize',12);
    plot(nmult_btfptj,err_btfptj,'o--','linewidth',2,'markersize',10);
    text(nmult_btfptj(idx_show_mult_btfptj)+50,...
        err_btfptj(idx_show_mult_btfptj)+0.8,...
        nlayers_btfptj(idx_show_mult_btfptj),'fontsize',12);
    
    axis([0,xrange_mult,-2,1.2*max(err_tj)]);
    xlabel('Number of multiplications per GFT');
    ylabel('Empirical averge error \epsilon');
    set(gca,'fontsize',18);
    legend('Matrix GFT', 'Haar-matrix-GFT', 'PTJ-GFT', 'Haar-PTJ-GFT');

    out_png=['figures/multvserr_', files{i}, '.png'];
    hgexport(gcf,out_png,hgexport('factorystyle'),'format','png');
    out_eps=['figures/multvserr_', files{i}, '.eps'];
    hgexport(gcf,out_eps,hgexport('factorystyle'),'format','eps');
    
    %% plot runtime vs empirical average error
    figure;
    plot(rt_mat/n_input,0,'*','linewidth',4,'markersize',20); 
    hold on; grid on;
    plot(rt_btf/n_input,0,'s','linewidth',4,'markersize',20);
    plot(rt_ptj/n_input,err_ptj,'+:','linewidth',2,'markersize',10);
    text(rt_ptj(idx_show_ptj)/n_input+1e-7,err_ptj(idx_show_ptj)+0.8,...
        nlayers_ptj(idx_show_ptj),'fontsize',12);
    plot(rt_btfptj/n_input,err_btfptj,'o--','linewidth',2,'markersize',10);
    text(rt_btfptj(idx_show_btfptj)/n_input+1e-7,err_btfptj(idx_show_btfptj)+0.8,...
        nlayers_btfptj(idx_show_btfptj),'fontsize',12);
    
    axis([0,xrange_rt,-2,1.2*max(err_tj)]);
    xlabel('Averge runtime per GFT (sec)');
    ylabel('Empirical averge error \epsilon');
    set(gca,'fontsize',18);
    legend('Matrix GFT', 'Haar-matrix-GFT', 'PTJ-GFT', 'Haar-PTJ-GFT');

    out_png=['figures/rtvserr_', files{i}, '.png'];
    hgexport(gcf,out_png,hgexport('factorystyle'),'format','png');
    out_eps=['figures/rtvserr_', files{i}, '.eps'];
    hgexport(gcf,out_eps,hgexport('factorystyle'),'format','eps');
    
    %% plot runtime vs sign-normalized RE
    switch files{i}
        case 'bd8x8'
            res_ptj=re_b_ptj_per_layer(1:5:end);
            res_btfptj=re_b_btfptj_per_layer(1:5:end);
        case 'z8x8'
            res_ptj=re_z_ptj_per_layer(1:5:end);
            res_btfptj=re_z_btfptj_per_layer(1:5:end);
    end
    figure;
    plot(rt_mat/n_input,0,'*','linewidth',4,'markersize',20); 
    hold on; grid on;
    plot(rt_btf/n_input,0,'s','linewidth',4,'markersize',20);
    plot(rt_ptj/n_input,res_ptj,'+:','linewidth',2,'markersize',10);
    text(rt_ptj(idx_show_ptj)/n_input+1e-7,res_ptj(idx_show_ptj)+4e-2,...
        nlayers_ptj(idx_show_ptj),'fontsize',12);
    plot(rt_btfptj/n_input,res_btfptj,'o--','linewidth',2,'markersize',10);
    text(rt_btfptj(idx_show_btfptj)/n_input+1e-7,res_btfptj(idx_show_btfptj)+4e-2,...
        nlayers_btfptj(idx_show_btfptj),'fontsize',12);
    
    axis([0,xrange_rt,-0.1,1.2*max(res_ptj)]);
    xlabel('Averge runtime per GFT (sec)');
    ylabel('Sign-normalized RE \delta');
    set(gca,'fontsize',18);
    legend('Matrix GFT', 'Haar-matrix-GFT', 'PTJ-GFT', ...
        'Haar-PTJ-GFT', 'location', 'northeast');

    out_png=['figures/rtvsre_', files{i}, '.png'];
    hgexport(gcf,out_png,hgexport('factorystyle'),'format','png');
    out_eps=['figures/rtvsre_', files{i}, '.eps'];
    hgexport(gcf,out_eps,hgexport('factorystyle'),'format','eps');
end
