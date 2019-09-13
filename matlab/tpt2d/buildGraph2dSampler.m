function [DB,T,resets] = buildGraph2dSampler(N)
    %create the graph and other properties by starting in a worm state,
    %having a strong backbone and weak other bonds, run BD sim for a long
    %time, keep track of states, connections, 
    
    %set parameters
    r = 40; %Range parameter to pair potential
    beta = 1; %inverse temp
    method = 1; % SDE solver. 1 -> EM. 2-> RK. 
    DT = 0.01;   %time increment after to check for bond formation
    Kh = 1850;    %large sticky parameter, for backbone
    Kl = 0.1;     %small sticky parameter, for other bonds
    state = 1;    %index of worm state 1
    timer = 0;    %keep track of time since last transition
    samples = 100000; %number of samples to run for
    resets = [];  %record times BD sim breaks
    
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
    El = stickyNewton(8,r,Kl);
    [X0, types, P, E] = setupBD(N, Eh, El);
    
    %do simulation samples number of times
    counter = 0;
    while counter < samples
        %run BD for DT seconds
        [~,Xf] = solveSDE(X0, N, DT, r, E, beta, P, method);
        
        %check if state is new or not
        [state, DB, T, timer, Xf, resets] = checkState(Xf, state, DB, T, timer, DT, resets);
        
        %update starting state
        X0 = Xf';
        counter = counter + 1
%         cPlot(Xf,types,P); 
%         pause()
    end
     
end

function [c, types, P, E] = setupBD(N, Eh, El)
    %create initial configuration and parameters to bd simulation
    %create worm
    c = zeros(2*N,1);
    c(1:2:2*N)=(1:N);
    %c(1:2:2*N) = (N:-1:1);
    
    %create types
    types = ones(N,1);    %all particles are the same type
    
    %create interactions
    P = ones(N);          %all particles bond
    
    %create energy matrix - strong backbone
    E = ones(N)*El; S = diag(diag(E,-1),-1); S = S + S';
    E = E-S; E = E+Eh/El*S;
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

function [new_state, DB, T, timer, X, resets] = checkState(X, prev_state, DB, T, timer, DT, resets)
    %check if the current state is changed from the previous.
    %if so, check if seen before and update quantities
    
    %first, update timer by DT
    timer = timer + DT; %gives current time
    
    %check if new state is connected
    oldA = DB{prev_state,1}; [newA, bonds] = getAdj(X); tol = 1e-5;
    N = length(diag(newA));
    if sum(diag(newA,1)) ~= N-1
        %core bonds broken, reset to starting config
        X = zeros(2*N,1)';
        X(1:2:2*N)=(1:N);
        timer = 0; new_state = 1;
        resets = [resets prev_state];
        return
    end
    
    %get old and new adj matrices, compare
    if sum(sum(abs(oldA-newA))) < tol
        %matrices are the same
        DB{prev_state,2} = DB{prev_state,2} + 1;
        new_state = prev_state;
        return
    else
        %matrices are not the same. Check if seen before by bonds
        [L,~] = size(DB);
        for i = 1:L
            if DB{i,3} == bonds
                %same number of bonds, compare matrices
                comp = DB{i,1}; 
                if sum(sum(abs(comp-newA))) < tol
                    %matrices are the same
                    new_state = i;
                    DB{new_state,2} = DB{new_state,2} + 1;
                    [sizea,~] = size(DB{new_state, 5});
                    if sizea < 100
                        DB{new_state, 5} = [DB{new_state, 5}; X];
                    end
                    %check if connection exists
                    if any(DB{new_state,4} == prev_state) == 0
                        %not connected, make 2 way connection
                        DB{new_state,4} = [DB{new_state,4} prev_state];
                        DB{prev_state,4} = [DB{prev_state,4} new_state];
                        T{prev_state,new_state}=[]; T{new_state, prev_state} = [];
                    end
                    T{prev_state,new_state} = [T{prev_state,new_state}, timer];
                    timer = 0;
                    return
                end
            end
        end
        %if this point is reached, newA is not in DB. Add it
        new_state = L+1;
        DB{new_state,5} = X;
        DB{new_state,1} = newA; DB{new_state,2} = 1; DB{new_state,3} = bonds;
        DB{new_state,4} = []; DB{new_state,4} = [DB{new_state,4} prev_state];
        DB{prev_state,4} = [DB{prev_state,4} new_state];
        T{prev_state,new_state}=[]; T{new_state, prev_state} = [];
        T{prev_state,new_state} = [T{prev_state,new_state}, timer];
        timer = 0;
        return
    end
end
    
    
    
    