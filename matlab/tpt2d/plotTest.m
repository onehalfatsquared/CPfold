N = 7;

%import the graph structure with all states
%Check if database exists for this value of N
filename = strcat('N', num2str(N), 'DB.mat');
if exist(filename, 'file')
    %file exists, load it
    load(filename);
else
    error('Requested database does not exist. Initialize it first');
end

cs = DB{1,5};