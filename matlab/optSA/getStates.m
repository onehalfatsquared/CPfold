function getStates(N, first)
    %This function runs long BD simulation starting in worm configuration.
    %Keeps track of clusters visited and checks against a database. Adds
    %new cluster to database (cell struct) if not there and increments a count if it
    %exists. Result simultaneously estimates partition fn of all clusters. 
    
    %Check if database exists for this value of N
    filename = strcat('N', num2str(N), 'allStates.mat'); 
    if exist(filename, 'file')
        %file exists, load it
        tempL = load(filename); %temporarily load struct 
        q = fieldnames(tempL);  %get field name in struct, only 1 field
        DB = getfield(tempL, {1}, q{1}); %get database
        clear tempL q; 
    else
        %file does not exist, create database
        DB = {}; 
    end
    
    %setup the BD sim
    [X0, types] = setupBD(N, first); %gives worm ic and types of particles
    K = [1,1,1]; %uniform sticky parameters
    P = ones(N); %interaction matrix - everything bonds here P=1. 
    rho = 40; %range parameter
    beta = 2; %inverse temp
    Tf = 1;   %final time, what to use?
    
    %get corresponding energy in morse potential, create pair energy matrix
    Ek = [stickyNewton(8,rho,K(1)), stickyNewton(8,rho,K(2)), stickyNewton(8,rho,K(3))];
    E = createE(types, Ek); 
    
    %run until state change. How often to check? Changes a lot at beginning
    %dont store whole trajectory, just current state and last adj matrix
    method = 1; 
    [t,Xf] = solveSDE(X0, N, Tf, rho, E, beta, P, method)
    cPlot(Xf', types, P); 
    
    %need adj matrix code
    
    
    %check if found before
    
    
end

function [c, types] = setupBD(N, first)
    %generate the initial condition for BD simulation - worm
    %types go ABAB or BABA until end. first specifies first type
    
    %create the worm, each particle dist 1 apart on x axis
    c = zeros(3*N,1);
    c(1:3:3*N) = (1:N);
    
    %create types vector
    reps = floor(N/2);  %number of repitions of AB/BA
    if first == 1
        types = repmat([1,2],1,reps);
    elseif first == 2
        types = repmat([2,1],1,reps);
    else
        error('Input for first must be 1 or 2')
    end
    %check if N is odd. if yes, then add first to end of types
    if rem(N,2) == 1
        types(N) = first;
    end
end

function [t,Xf] = solveSDE(X0,N,T,rho,E,beta,P,method)
    %solve SDE with EM. Output X(T)
    % method 1 -> EM. method 2 -> RK scheme. 
    if method == 1
        k = 1e-6; Nt = round(T/k); %time step, num time steps
        [t,Xf] = EM(X0,N,k,Nt,rho,E,beta,Nt,P);
    elseif method == 2 %Not sure if working correctly. 
        k = 1e-8; Nt = round(T/k); %time step, num time steps
        [t,Xf] = RK(X0,N,k,Nt,rho,E,beta,Nt,P); 
    end
end
    

function DB = findAdj(A, DB)
    %checks if adjacency matrix A is present in the database DB.
    %increment count if yes, add if no. 
    tol = 1e-6; 
    for i = 1:length(DB)
        B = DB{i,1}; 
        if sum(sum(abs(A-B))) < tol
            %matrices are same, increment count
            DB{i,2} = DB{i,2} + 1;
            fprintf('Matrix %d found again. Count: %d', i, DB{i,2});
            break
        end
    end
    %matrix does not exist yet, add
    fprintf('New state found. Count: %d', i+1);
    DB{i+1,1} = A; 
    DB{i+1,2} = 1;
end