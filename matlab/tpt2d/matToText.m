N = 6;

%import the graph structure with all states
%Check if database exists for this value of N
filename = strcat('N', num2str(N), 'DB.mat');
if exist(filename, 'file')
    %file exists, load it
    load(filename);
else
    error('Requested database does not exist. Initialize it first');
end

outfile = fopen(strcat('N', num2str(N), 'DB.txt'),'w');
fprintf(outfile, strcat(num2str(N),'\n'));
fprintf(outfile, strcat(num2str(length(DB)),'\n'));
for i = 1:length(DB)
    A = DB{i,1};
    for q = A(:)
        fprintf(outfile, "%d ", q);
    end
    fprintf(outfile, "%d ", DB{i,2});
    fprintf(outfile, "%d ", DB{i,3});
    coord = DB{i,5}; [m,~] = size(coord); coord=coord';
    fprintf(outfile, "%d ", m);
    for q = coord(:)
        fprintf(outfile, "%f ", q);
    end
    fprintf(outfile, "N\n");
end
fclose(outfile);
