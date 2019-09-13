function [t,T] = RK(X0,N,k,Nt,rho,E,beta,SS,P)
    %solve SDE with RK. Subsample every SS timesteps. Output trajectory. 
     
    X = X0; sample = 1; 
    %if subsampling, set IC as first sample.
    if SS~= Nt
        L = round(Nt/SS)+1;  %Number of sample points
        T = zeros(L,2*N); t = zeros(1,L); %initialize storage
        T(1,:) = X0'; t(1) = 0; %set first sample to IC
        sample = 2;  
    end
    
    %time-stepping using RK scheme
    for i=1:Nt
        %compute arguments for predictor step 
        xi = randn(2*N,1)*sqrt(1/beta);  %Normal rv with 1/beta var
        Xtilde = X + sqrt(k/2)*xi;       %first eval point
        particles = c2p(Xtilde);         %particle array
        
        %predictor step, evaluate Gk, get Xstar
        g0 = mGrad(particles, rho, E, N, P);       %force at Xtilde
        x1 = Xtilde - 2/3*k*g0;                    %compute step
        g1 = mGrad(c2p(x1), rho, E, N, P);         %force at x1
        Gk = -1/4*g0 -3/4*g1;                      %RK formula
        Xstar = 2*Xtilde - X +k*Gk;                %proposed move       
        
        %corrector step, compute acc prob, update
        eta = xi + sqrt(2*k)*Gk;
        U0 = mEval(c2p(X), rho, E, N, P);
        U1 = mEval(c2p(Xstar), rho, E, N, P);
        alpha = min(1,exp(-beta*((eta'*eta-xi'*xi)/2+U1-U0))) %acc prob
        gamma = (rand<alpha);                      %bernouli rv
        %gamma = 1;                                 %???? if not, alpha = 0..
        X = gamma*Xstar + (1-gamma)*X;             %update positions
        
        %subsample
        if mod(i,SS)==0
            T(sample,:) = X'; t(sample) = (i-1)*k;
            sample = sample+1;
        end
    end
end
