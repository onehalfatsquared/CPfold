function [tau] = estimateMFPT(N, state)
    %estimate the mean first passage time from state i to j, where state j
    %has more bonds than state i (time for diffusion to bring particles in
    %range). Equilibrates trajectory in state i, resets when state j is
    %reached. 
    
    %set parameters
    r = 40; %Range parameter to pair potential
    beta = 1; %inverse temp
    method = 1; % SDE solver. 1 -> EM. 2-> RK. 
    DT = 0.01;   %time increment after to check for bond formation
    Kh = 1850;    %large sticky parameter, for backbone
    timer = 0;    %keep track of time since last transition
    samples = 25; %number of samples to run for
    eq = 20;     %number of equilibration steps
    tau = cell(700); % to be replaced by a loaded cell
    
    %import the graph structure and transition times
    %Check if database exists for this value of N
    filename = strcat('N', num2str(N), 'DB.mat'); 
    if exist(filename, 'file')
        %file exists, load it
        load(filename);
    else
        error('Requested database does not exist. Initialize it first');
    end
  
    %setup simulation
    Eh = stickyNewton(8,r,Kh);
    [start, types, P, E] = setupBD(N, Eh);
    
    %do simulation samples number of times
    counter = 0; 
    while counter < samples
        %equilibrate the trajectory
        eq_count = 0; new_state = 0; reflect = 0; X0 = start;
        while eq_count < eq
            %run BD for DT seconds
            [~,Xf] = solveSDE(X0, N, DT, r, E, beta, P, method);
            %check if state changed
            [new_state, timer, reflect, reset] = checkState(Xf, state, DB, timer, DT);
            if reflect == 0 && reset == 0
                X0 = Xf';
            end
            eq_count = eq_count + 1;
        end
        %trajectory equilibrated, do test
        timer = 0 ; reflect = 0;
        while reflect == 0
            %run BD for DT seconds
            [~,Xf] = solveSDE(X0, N, DT, r, E, beta, P, method);
            %check if state changed
            [new_state, timer, reflect, reset] = checkState(Xf, state, DB, timer, DT);
            if reflect == 0 && reset == 0
                X0 = Xf';
            elseif reflect == 1
                tau{state, new_state} = [tau{state,new_state} timer];
            end
        end
        counter = counter + 1
    end
end

function [c, types, P, E] = setupBD(N, Eh)
    %create initial configuration and parameters to bd simulation
    %create worm
    c = zeros(2*N,1);
    c(1:2:2*N)=(1:N);
    
    %create types
    types = ones(N,1);    %all particles are the same type
    
    %create interactions
    P = ones(N);          %all particles bond
    
    %create energy matrix - strong bonds
    E = ones(N)*Eh;
end

function [t,T] = solveSDE(X0,N,T,rho,E,beta,P,method)
    %solve SDE with EM. Subsample every SS timesteps. Output trajectory. 
    % method 1 -> EM. method 2 -> RK scheme. 
    if method == 1
        k = 2.5e-6; Nt = round(T/k); %time step, num time steps
        SS = Nt; %Sub-sampling interval
        [t,T] = EM(X0,N,k,Nt,rho,E,beta,SS,P); 
    elseif method == 2 %Not sure if working correctly. 
        k = 1e-8; Nt = round(T/k); %time step, num time steps
        SS = Nt; %Sub-sampling interval
        [t,T] = RK(X0,N,k,Nt,rho,E,beta,SS,P); 
    end
end

function [A, b] = getAdj(clust)
    %get the adjacency matrix of a cluster and number of bonds present
    
    %create particle array, get distance matrix
    particles = c2p(clust); L=length(particles); D = zeros(L);
    for i=1:L
        for j=i+1:L
            D(i,j)=euDist(particles,i,j);
        end
    end
    
    %Initialize adjacency matrix and compute
    tol = 1e-5; A = zeros(L); b = 0;
    for i=1:L
        for j=i+1:L
            if D(i,j) < 1.1+tol % range of potential for rho=38
                b = b+1; 
                A(i,j)=1; A(j,i)=1;
            end
        end
    end
end

function [new_state, timer, reflect, reset] = checkState(X, prev_state, DB, timer, DT)
    %check if the current state is changed from the previous.
    
    %first, update timer by DT
    timer = timer + DT; %gives current time
    reflect = 0; reset = 0;
    
    %check if new state is connected
    oldA = DB{prev_state,1}; [newA, bonds] = getAdj(X); tol = 1e-5;
    N = length(diag(newA));
    if sum(diag(newA,1)) ~= N-1
        %core bonds broken, reflect to starting config
        new_state = 0;
        reset = 1 
        return
    end
    
    %get old and new adj matrices, compare
    if sum(sum(abs(oldA-newA))) < tol
        %matrices are the same
        new_state = prev_state;
        return
    else
        %matrices are not the same. Find new state
        [L,~] = size(DB);
        for i = 1:L
            if DB{i,3} == bonds
                %same number of bonds, compare matrices
                comp = DB{i,1}; 
                if sum(sum(abs(comp-newA))) < tol
                    %matrices are the same
                    new_state = i;
                    reflect = 1;
                    return
                end
            end
        end
    end
end
    
    
    
    