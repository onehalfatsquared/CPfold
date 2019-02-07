function stickyE=stickyNewton(E,rho,kappa0)
    %Uses root finding algo to find energy correspodning to given value of
    %kappa. 
    
    x0=E; cutoff=100; tol=1e-4; %Parameters
    
    if x0>0.7 || x0<0.3
        %Perform newtons method
        for step=1:cutoff
            [a,b]=stickyF(x0, rho, kappa0); %f and f'
            x1=x0-a/b;   %Newton step
            if abs(x1-x0)<tol
                stickyE=x1;
                return
            end
            x0=x1;
        end
        stickyE=x1;
    else
        %perform biscetion method
        a=0.1; b=E; %endpoints of bisection interval
        fprintf('Bisection\n')
        for step=1:cutoff
            c=(a+b)/2; f=stickyF(c,rho,kappa0); %midpoint and fn val
            if abs(f)<tol
                stickyE=c;
                return
            end
            if f*stickyF(a,rho,kappa0)>0
                a=c;
            else
                b=c;
            end
        end
        stickyE=c;
    end
end