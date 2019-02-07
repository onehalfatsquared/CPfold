function DB = buildGraph2d(N)
    %This function builds a database of adjacency matrices for transition
    %states of N spheres in 2d. Starts with a rigid state and breaks bonds
    %one at a time to get new states. Stores a state index, and adjacency
    %matrix, list of connected states, and number of bonds. (eventually
    %coordinates?) 
    
    %Check if database exists for this value of N
    filename = strcat('N', num2str(N), '.mat'); 
    if exist(filename, 'file')
        %file exists, load it
        tempL = load(filename); %temporarily load struct 
        q = fieldnames(tempL);  %get field name in struct, only 1 field
        DB = getfield(tempL, {1}, q{1}); %get database
        clear tempL q; 
    else
        %file does not exist, create database
        DB = {}; 
        
        %add the starting rigid states
        if N == 7
            %state 1
            DB{1,1} = [0,1,0,0,0,1,1;
                       1,0,1,0,0,0,1;
                       0,1,0,1,0,0,1;
                       0,0,1,0,1,0,1;
                       0,0,0,1,0,1,1;
                       1,0,0,0,1,0,1;
                       1,1,1,1,1,1,0];
           DB{1,2} = [];
           DB{1,3} = 12;
           
           %state 2
           DB{2,1} = [0,1,1,0,0,0,0;
                      1,0,1,1,0,0,0;
                      1,1,0,1,0,0,0;
                      0,1,1,0,1,1,0;
                      0,0,1,1,0,1,1;
                      0,0,0,1,1,0,1;
                      0,0,0,0,1,1,0];
           DB{2,2} = [];
           DB{2,3} = 11;
           
           %state 3
           DB{3,1} = [0,1,1,0,0,0,0;
                      1,0,1,1,0,0,0;
                      1,1,0,1,1,0,0;
                      0,1,1,0,1,1,1;
                      0,0,1,1,0,1,0;
                      0,0,0,1,1,0,1;
                      0,0,0,1,0,1,0];
           DB{3,2} = [];
           DB{3,3} = 11;
                      
           %state 4
           DB{4,1} = [0,1,1,0,0,0,0;
                      1,0,1,0,0,1,1;
                      1,1,0,1,0,1,0;
                      0,0,1,0,1,1,0;
                      0,0,0,1,0,1,0;
                      0,1,1,1,1,0,1;
                      0,1,0,0,0,1,0];
           DB{4,2} = [];
           DB{4,3} = 11;
           
           %state 5
           
        end
    end
    
    %get the max number of bonds in any state so far
    [L, ~] = size(DB);
    M = 0;
    for i = 1:L
        m = DB{i,3};
        if m >= M
            M = m;
        end
    end
    
    %make vector of allowed number of bonds, N-1 to M, backwards
    num_bonds = M:-1:N-1;
    
    %create indices of 1-7 chain that cant be broken -> superdiagonal,
    %column indexed
    chain = zeros(1,N-1);
    for i = 1:N-1
        chain(i) = i*(N+1);
    end
    
    %Loop over levels (number of broken bonds, 0 to N-2)
    for levels = num_bonds
        [L,~] = size(DB) %get length
        for i = 1:L
            %until the end of the existing database, check if 
            %level is = DB{i,3}. While it is, take the matrices 
            %and break bonds
            if DB{i,3} == levels
                A = DB{i,1}; %adjacency matrix
                Au = triu(A); %upper triangular part
                
                %get index of bonds - exclude the chain
                ind = find(Au>0); ind = setdiff(ind,chain);
                
                %break the bonds one by one and check if 
                %permissible. if yes, check if present. 
                for j = 1:length(ind)
                    B = Au;
                    B(ind(j)) = 0; B = B + B';
                    c = checkP(B, N); 
                    if c == 1
                        [DB] = findAdj(i, B, DB, L, levels-1);
                    end
                end
            end
        end
    end 
end
    

function [DB] = findAdj(parent, B, DB, oldL, levels)
    %checks if adjacency matrix B (broken) is present in the database DB.
    %If not add it, and give a connection to A (parent). If yes, check if 
    %connected to A and connect if not
    
    tol = 1e-6;  %tolerance for diff in matrices elementwise
    [newL,~] = size(DB); %number of matrices at current level
    diff = newL - oldL; %number of matrices to compare to
    test = oldL;   %for the case of diff = 0
    for i = 1:diff
        test = i+oldL; %index of matrices to compare to
        comp = DB{test,1}; %matrix to compare to
        if sum(sum(abs(B-comp))) < tol
            %matrices are same, check connection
            if any(DB{parent,2} == test) == 0
                %not connected, make 2 way connection
                DB{parent,2} = [DB{parent,2} test];
                DB{test,2} = [DB{test,2} parent];
            end
            return
        end
    end
    %matrix does not exist yet, add, make connection
    DB{test+1,1} = B; 
    DB{test+1,2} = [DB{test+1,2} parent];
    DB{parent,2} = [DB{parent,2} test+1];
    DB{test+1,3} = levels; 
end

function [c] = checkP(A, N)
    %perform dfs to check graph is connected. return 1 if yes, 0 if no
    
    %create a vector of visited states, and states to visit. start at 1
    s = 1; visited = zeros(1,N); toVisit=s; 
    
    %loop while states to visit in non empty
    while isempty(toVisit) == 0
        visited(toVisit(1)) = 1;
        C = A(toVisit(1),:);
        children = find(C>0);
        toVisit(1)=[];
        for child = children
            if visited(child) == 0
                toVisit = [toVisit child];
            end
        end
    end
    
    %if all states have been visited, graph is connected
    if visited == ones(1,N)
        c = 1;
    else 
        c = 0;
    end
end

function [b] = checkPerm(A, B, N)
    %function to check if the adjacency matrices are permutations of each
    %other. Starts by ruling out matrices with different degree
    %distributions, then checks all perms. 
    
    %compute the sorted list of degrees
    d1 = sort(sum(A)); d2 = sort(sum(B)); 
    tol = 1e-4; 
    
    %check if degree lists are same, if not return 0
    b=0; 
    if d1 ~= d2
        return
    else
        %construct all permutations and check for isomorphism
        allPerms = perms(1:N); 
        I = eye(N);
        for i = 1:length(allPerms)
            P = I(allPerms(i,:),:);
            C = P*A*P';
            S = sum(sum(abs(C-B)));
            if S < tol
                b = 1;
                return
            end
        end
    end   
end