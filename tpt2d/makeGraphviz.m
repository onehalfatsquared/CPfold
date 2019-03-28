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

sum = 0;
%get total number of trials
for i = 1:length(DB)
    sum = sum+DB{i,2};
end

%open the file to write to
d = 1; prev = 1; width = 1;
filename=strcat('N',num2str(N),'tptGraph.txt'); %File name 
fileID=fopen(filename,'w');                     %Open file
fprintf(fileID, 'graph tpt_graph%d  {\n nodesep=0.1;\n',N); %header to start graph environment

%make the initial cluster seperately, no edges
fprintf(fileID,'"c%d" [label="%d"  ,  shape=circle  , width=%f, regular=1,style=filled,fillcolor=white] ;\n', d, d, 4);      %Write cluster
fprintf(fileID,'{rank = same; "c1";}');
%fclose(fileID);


%make the connections starting at clusters with 7 bonds

for b = 7:12
    for i = 1:length(DB)
        first = 0; added=[];
        %check if state has b bonds
        if DB{i,3} == b
            %get the list of connections, make connections to existing
            %nodes if necc
            conn = DB{i,4};
            for j = conn
                if DB{j,3} == b-1 
                    %first time through make the node
                    if first == 0
                        %filename=strcat('N',num2str(N),'tptGraph.txt'); %File name 
                        %fileID=fopen(filename,'w');                     %Open file
                        freq = DB{i,2}; width = freqToWidth(freq, sum)
                        added=[added i];
                        fprintf(fileID,'"c%d" [label="%d"  ,  shape=circle  , width=%f, regular=1,style=filled,fillcolor=white] ;\n', i, i, width);%Write cluster
                        fprintf(fileID,'"c%d" -- "c%d" ;\n', j, i);
                        first = 1; prev = i;
                        %fclose(fileID);
                    else
                        %after just make edges
                        fprintf(fileID,'"c%d" -- "c%d" ;\n', j, i);
                    end
                end
            end
        end
    end
    fprintf(fileID,'{rank = same; ');
    for a = 1:length(added)
        fprintf(fileID,'"c%d";', a);
    end
    fprintf(fileID,'}\n');
end

filename=strcat('N',num2str(N),'tptGraph.txt'); %File name 
fileID=fopen(filename,'a+');                     %Open file
fprintf(fileID, '}');                           %end the graphviz file with }
fclose(fileID);                                 %Close file




function w = freqToWidth(freq, sum)
    if freq/sum < 1/800
        w = 1;
    else
        w = 1+freq/sum*32;
    end
end