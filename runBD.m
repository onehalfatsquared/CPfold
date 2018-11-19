%Run BD simulation - dXt = -grad(U(Xt))+sqrt(2/beta)dWt
clear
%Cool tests: N=6: 11,210,1007
%            N=9: 12, 33, 
%            N=10: 12, 

%set parameters
rng(33) %set random seed   
N = 9; %Number of particles
r = 6; %Range parameter to pair potential
E = 5; %Well-depth to pair potential
beta = 10; %inverse temp
T = 6; %final time
k = 2e-4; Nt = T/k; %time step, num time steps
SS = 300; %Sub-sampling interval

%setup simulation
[X0,types] = setIC(N,'rand'); %set initial configuration
P = createPI(types);          %create interaction matrix

%do simulation
[t,T] = solveSDE(X0,N,k,Nt,r,E,beta,SS,P); %get a trajectory
Xf = T(end,:)'; %get final state

%make plots
cPlot(Xf,types,P); 
makeMovie(T,types,P); 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c,types] = setIC(N,type)
    %Set initial condition. particles in a straight line - worm
    c = zeros(3*N,1);
    c(1:3:3*N)=(1:N);
    %Set distribution of particle types. 3 kinds of particles.
    if type == 'same'
        types = ones(N,1);
    elseif type == 'rand'
        types = randi(3,N,1);
    end
end


function [t,T] = solveSDE(X0,N,k,Nt,rho,E,beta,SS,P)
    %solve SDE with EM. Subsample every SS timesteps. Output trajectory. 
    X = X0; sample = 1; 
    %if subsampling, set IC as first sample.
    if SS~= Nt
        T(1,:) = X0'; t(1) = 0;
        sample = 2;
    end
    for i=1:Nt
        particles = c2p(X); %particle array
        a = -morseGrad(particles,rho,E,N,P)*k; %det part
        b = randn(3*N,1)*sqrt(2/beta*k); %stoch part
        X = X + a + b; %EM step
        %subsample
        if mod(i,SS)==0
            T(sample,:) = X'; t(sample) = (i-1)*k;
            sample = sample+1;
        end
    end
end

function makeMovie(T,types,P)
    %make movie of evolution of cluster
    [N,~]=size(T);
    for i=1:N
        X = T(i,:)'; %cluster at time t_i
        cPlot(X,types,P);  
    end
end
    




