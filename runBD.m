%Run BD simulation - dXt = -grad(U(Xt))+sqrt(2/beta)dWt
clear
%Cool tests, rng3:  N=6: 11,210,1007
%                   N=9: 12, 33, 
%                   N=10: 12, 

try mex mGrad.cpp; catch,  end

%set parameters
rng(33) %set random seed   
N = 9; %Number of particles
r = 6; %Range parameter to pair potential
K = [20,10,20]; %vector of sticky parameters. (1->1, 1<->2, 2->2). 
beta = 10; %inverse temp
T = 3; %final time
k = 2e-4; Nt = T/k; %time step, num time steps
SS = 300; %Sub-sampling interval
method = 1; % SDE solver. 1 -> EM. 2-> RK. 

%setup simulation
[X0,types] = setIC(N,'rng2'); %set initial configuration
P = createPI(types);          %create interaction matrix
%Do newton solve for energies, then compute matrix

%do simulation
[t,T] = solveSDE(X0,N,k,Nt,r,E,beta,SS,P,method); %get a trajectory
Xf = T(end,:)'; %get final state

%make plots
cPlot(Xf,types,P); 
M = makeMovie(T,types,P); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c,types] = setIC(N,type)
    %Set initial condition. particles in a straight line - worm
    c = zeros(3*N,1);
    c(1:3:3*N)=(1:N);
    %Set distribution of particle types. 
    if type == 'same'
        types = ones(N,1);
    elseif type == 'rng3'
        types = randi(3,N,1);
    elseif type == 'rng2'
        types = randi(2,N,1);
    end
end

function [t,T] = solveSDE(X0,N,k,Nt,rho,E,beta,SS,P,method)
    %solve SDE with EM. Subsample every SS timesteps. Output trajectory. 
    % method 1 -> EM. method 2 -> RK scheme. 
    if method == 1
        [t,T] = EM(X0,N,k,Nt,rho,E,beta,SS,P);
    elseif method == 2
        [t,T] = RK(X0,N,k,Nt,rho,E,beta,SS,P); 
    end
end

function M = makeMovie(T,types,P)
    %make movie of evolution of cluster
    [N,~]=size(T);
    M(N) = struct('cdata',[],'colormap',[]);
    for i=1:N
        X = T(i,:)'; %cluster at time t_i
        cPlot(X,types,P);  
        M(i) = getframe(gcf); 
    end
end
    




