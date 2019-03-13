%make graphviz representation of the state data. COnnections between
%adjacnecy matrices with transitions, size of node preoprtional to
%frequency. 

clear;

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


%make the initial cluster seperately, no edges



%make the connections starting at clusters with 7 bonds

for b = 7:12 
    for i = 1:length(DB)
        %check if state has b bonds
        if DB{i,3} == b
            %get the list of connections, make connections to existing
            %nodes if necc
            conn = DB{i,4};
            for j = conn
                if DB{j,3} == b-1 
                    %first time through make the node
                    
                    %after just make edges
                    
                end
            end
        end
    end
end
