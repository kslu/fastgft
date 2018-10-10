function print_c_gftmatrix(U,txname,filename)
n=size(U,1);
fid=fopen(filename, 'w+');

fwrite(fid, ['static const double ', txname, '[', sprintf('%d',n^2), ...
    sprintf('] = {\n')], 'char');

% write 5 values per line
num_per_line=5;
for j=1:ceil(n^2/num_per_line)
    
    for k=(j-1)*num_per_line+1:min(n^2,j*num_per_line)
        fwrite(fid, sprintf('%.10f', U(k)));
        if k~=n^2
            fwrite(fid, ', ');
        end
    end
    fwrite(fid, sprintf('\n'));

end

fwrite(fid, sprintf('};\n\n'));