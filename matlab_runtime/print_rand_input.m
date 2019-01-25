function print_rand_input(n,m,filename)
% print_rand_input(n,m,filename)
% n: transform length
% m: number of blocks
fid=fopen(filename, 'w+');
fwrite(fid, sprintf('%d \n', m));

X=rand(m,n);

for i=1:m
    for j=1:n
        fwrite(fid, sprintf('%.10f ', X(i,j)));
    end
    fwrite(fid, sprintf('\n'));
end