function [t,T] = EM(X0,N,k,Nt,rho,E,beta,SS,P)
    %solve SDE with EM. Subsample every SS timesteps. Output trajectory. 
     
    X = X0; sample = 1; 
    %if subsampling, set IC as first sample.
    if SS~= Nt
        L = round(Nt/SS)+1;  %Number of sample points
        T = zeros(L,2*N); t = zeros(1,L); %initialize storage
        T(1,:) = X0'; t(1) = 0; %set first sample to IC
        sample = 2;  
    end
    
    %time-stepping using EM scheme
    for i=1:Nt
        particles = c2p(X); %particle array
        %a = -morseGrad(particles,rho,E,N,P)*k; %det part
        a = -mGrad(particles, rho, E, N, P)*k;
        
        b = randn(2*N,1)*sqrt(2/beta*k); %stoch part
        X = X + a + b; %EM step
        %subsample
        if mod(i,SS)==0
            T(sample,:) = X'; t(sample) = (i-1)*k;
            sample = sample+1;
        end
    end
end

