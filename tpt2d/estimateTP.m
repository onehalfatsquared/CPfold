function [out] = estimateTP(N, a)
    %estimate the transition probabilities for a given cluster forming bonds
    %through a monte carlo type method
    
    %set parameters
    r = 40; %Range parameter to pair potential
    e = 10;  %energy parameter
    beta = 1; %inverse temp
    method = 1; % SDE solver. 1 -> EM. 2-> RK. 
    M = 1000;   %number of monte carlo trials
    DT = 0.1;   %time increment after to check for bond formation
    
    %import the graph structure - may change to load and update
    G = buildGraph2d(N);
    
    %initialize the transition probability matrix - may change to load
    P = zeros(length(G));

    %setup simulation
    types = ones(N,1);    %all particles are the same type
    P = ones(N);          %create interaction matrix - all bonding
    E = e*ones(N);        %create energy matrix - all same
    
    
    %loop over the elements of DB backwards doing MC
    for i=0:length(G)-1
        %get all the relevant info
        coords = G{end-i,4};  %gets coordinates of the point
        NB = G{end-i,3};  %gets the number of bonds in starting state
        links = G{end-i,2}; %gets possibilities for current state to become
        
        %do stuff
        X0 = coords;
        [t,Xf] = solveSDE(X0, N, DT, r, E, beta, P, method);
        
        %do other stuff
        
    end
    
    
    
end




function [t,T] = solveSDE(X0,N,T,rho,E,beta,P,method)
    %solve SDE with EM. Subsample every SS timesteps. Output trajectory. 
    % method 1 -> EM. method 2 -> RK scheme. 
    if method == 1
        k = 1e-5; Nt = round(T/k); %time step, num time steps
        SS = Nt; %Sub-sampling interval
        [t,T] = EM(X0,N,k,Nt,rho,E,beta,SS,P);
    elseif method == 2 %Not sure if working correctly. 
        k = 1e-8; Nt = round(T/k); %time step, num time steps
        SS = Nt; %Sub-sampling interval
        [t,T] = RK(X0,N,k,Nt,rho,E,beta,SS,P); 
    end
end