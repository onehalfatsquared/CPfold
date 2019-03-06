%Run BD simulation - dXt = -grad(U(Xt))+sqrt(2/beta)dWt
clear

%Cool tests, rng3:  N=6: 11,210,1007
%                   N=9: 12, 33, 
%                   N=10: 12, 

%set parameters
%rng(33) %set random seed   
N = 7; %Number of particles
r = 40; %Range parameter to pair potential
K = [1,1,1]; %vector of sticky parameters. (1->1, 1<->2, 2->2). 
beta = 1; %inverse temp
T = 0.01; %final time
frames = 61;  %Number of frames in visualization
method = 1; % SDE solver. 1 -> EM. 2-> RK. 

%setup simulation
[X0,types] = setIC(N,'same'); %set initial configuration
P = ones(N);          %create interaction matrix
E = ones(N)*10; 

%do simulation
[t,T] = solveSDE(X0,N,T,r,E,beta,P,method,frames); %get a trajectory
Xf = T(end-1,:)'; %get final state

%make plots
cPlot(Xf,types,P); 
%M = makeMovie(t,T,types,P); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c,types] = setIC(N,type)
    %Set initial condition. particles in a straight line - worm
    c = zeros(2*N,1);
    c(1:2:2*N)=(1:N);
    %Set distribution of particle types. 
    if type == 'same'
        types = ones(N,1);
    elseif type == 'rng3'
        types = randi(3,N,1);
    elseif type == 'rng2'
        types = randi(2,N,1);
    end
end

function [t,T] = solveSDE(X0,N,T,rho,E,beta,P,method,frames)
    %solve SDE with EM. Subsample every SS timesteps. Output trajectory. 
    % method 1 -> EM. method 2 -> RK scheme. 
    if method == 1
        k = 5e-6; Nt = round(T/k); %time step, num time steps
        SS = round(Nt/frames); %Sub-sampling interval
        [t,T] = EM(X0,N,k,Nt,rho,E,beta,SS,P);
    elseif method == 2 %Not sure if working correctly. 
        k = 1e-8; Nt = round(T/k); %time step, num time steps
        SS = round(Nt/frames); %Sub-sampling interval
        [t,T] = RK(X0,N,k,Nt,rho,E,beta,SS,P); 
    end
end

function M = makeMovie(t,T,types,P)
    %make movie of evolution of cluster
    [N,~]=size(T);
    M(N) = struct('cdata',[],'colormap',[]);
    for i=1:N
        X = T(i,:)'; %cluster at time t_i
        cPlot(X,types,P);  
        title(strcat('t= ',num2str(t(i))),'fontsize',18);
        M(i) = getframe(gcf); 
    end
end
    




