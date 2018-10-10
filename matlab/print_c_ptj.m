function print_c_ptj(L,txname,coords_file,angle_file,idx_file,num_layers)
% print_c_ptj(L,txname,cfilename,afilename,err_thr)
if nargin<6
    num_layers=100;
end
if nargin<5 || isempty(idx_file)
    idx_file = [txname, '_idx.h'];
end
if nargin<4 || isempty(angle_file)
    angle_file=[txname, '_angles.h'];
end
if nargin<3 || isempty(coords_file)
    coords_file=[txname, '_coords.h'];
end

n=size(L,1);
[facs,~,~,~,~,idx_freq]=givens_ptj(L,num_layers*n/2,n/2);
[~,idx_orig]=sort(diag(L)');
idx_freq_c=[idx_orig; idx_freq]-1;  % minus 1 to match indexing in C

%% print c code
n_layers=length(facs);
n_per_line_coords=8;
n_per_line_double=2;
n_per_line_idx=15;
fcid=fopen(coords_file, 'w+');
faid=fopen(angle_file, 'w+');
fiid=fopen(idx_file, 'w+');

% print coordinates
fwrite(fcid, ['static const int ', txname, '_coords[', ...
    sprintf('%d * %d * 2] = {\n',n_layers,n/2)], 'char');

% print cos and sin
fwrite(faid, ['static const double ', txname, '_angles[', ...
    sprintf('%d * %d * 2] = {\n',n_layers,n/2)], 'char');

% print order of GFT coefficients
fwrite(fiid, ['static const int ', txname, '_idx[', ...
    sprintf('%d * %d] = {\n  ', num_layers+1, n)]);

for i=1:n_layers
    
    % parse the structure array
    cur_fac=facs{i};
    [i1,i2]=ind2sub([n,n],find(triu(facs{i},1)~=0));
    % include the identity part (nodes that are not paired)
    i_rest=setdiff(1:n,[i1;i2]);
    i1=[i1',i_rest(1:2:end)];
    i2=[i2',i_rest(2:2:end)];
    
    anglecos=full(cur_fac(sub2ind([n,n],i1,i1)));
    anglesin=full(cur_fac(sub2ind([n,n],i1,i2)));
    
    fwrite(fcid, '  ');
    fwrite(faid, '  ');
    for j=1:n/2
        if cur_fac(i1(j),i2(j))==0
            fwrite(fcid, '-1, -1, ');  % skip this Givens rotation
        else
            fwrite(fcid, sprintf('%d, %d, ',i1(j)-1,i2(j)-1));
        end
        if mod(j,n_per_line_coords)==0 && j~=n/2
            fwrite(fcid, sprintf('\n  '));
        end
        
        fwrite(faid, sprintf('%.10f, %.10f, ',anglecos(j),anglesin(j)));
        if mod(j,n_per_line_double)==0 && j~=n/2
            fwrite(faid, sprintf('\n  '));
        end
    end
    fwrite(fcid, sprintf('\n'));
    fwrite(faid, sprintf('\n'));
end

for i=1:size(idx_freq_c,1)
    for j=1:n
        fwrite(fiid, sprintf('%d, ', idx_freq_c(i,j)));
        if mod(j, n_per_line_idx)==0 && j~=n
            fwrite(fiid, sprintf('\n  '));
        end
    end
    
    fwrite(fiid, sprintf('\n  '));
end

fwrite(fcid, sprintf('};\n\n'));
fwrite(faid, sprintf('};\n\n'));
fwrite(fiid, sprintf('};\n\n'));