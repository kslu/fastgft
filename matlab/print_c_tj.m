function print_c_tj(L,txname,coord_file,angle_filename,idx_file,num_givens)
% not finished yet
if nargin<6
    num_givens=1000;
end
if nargin<5 || isempty(idx_file)
    idx_file=[txname, '_idx.h'];
end
if nargin<4 || isempty(angle_filename)
    angle_filename=[txname, '_angles.h'];
end
if nargin<3 || isempty(coord_file)
    coord_file=[txname, '_coords.h'];
end

n=size(L,1);
[fac_tj,~,~,~,~,idx_freq]=givens_tj(L,num_givens);
[~,idx_orig]=sort(diag(L)');
idx_freq_c=[idx_orig; idx_freq]-1;  % minus 1 to match indexing in C

% identify the "real" number of Givens rotations found
num_givens_found=num_givens;
for i=1:num_givens
    if nnz(triu(fac_tj{i},1))==0
        num_givens_found=i-1;
        break;
    end
end

%% print c code
n_per_line_coords=8;
n_per_line_double=2;
n_per_line_idx=15;
fcid=fopen(coord_file, 'w+');
faid=fopen(angle_filename, 'w+');
fiid=fopen(idx_file, 'w+');

% print coordinates
fwrite(fcid, ['static const int ', txname, '_coords[', ...
    sprintf('%d * 2 + 2] = {\n  ', num_givens_found)], 'char');

% print cos and sin
fwrite(faid, ['static const double ', txname, '_angles[', ...
    sprintf('%d * 2] = {\n  ', num_givens_found)], 'char');

% print order of GFT coefficients
fwrite(fiid, ['static const int ', txname, '_idx[', ...
    sprintf('%d * %d] = {\n  ', num_givens+1, n)]);

for i=1:num_givens_found
    
    % parse the structure array
    cur_fac=fac_tj{i};
    [i1,i2]=ind2sub([n,n],find(triu(cur_fac,1)~=0)); % should only have 1 element
    
    anglecos=full(cur_fac(i1,i1));
    anglesin=full(cur_fac(i1,i2));
    
    fwrite(fcid, sprintf('%d, %d, ',i1-1,i2-1));
    if mod(i,n_per_line_coords)==0
        fwrite(fcid, sprintf('\n  '));
    end

    fwrite(faid, sprintf('%.10f, %.10f, ',anglecos,anglesin));
    if mod(i,n_per_line_double)==0
        fwrite(faid, sprintf('\n  '));
    end
end

% write -1, -1 at the end for early stopping
fwrite(fcid, sprintf('-1, -1, \n'));

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