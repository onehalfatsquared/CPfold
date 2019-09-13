function [f,fprime] = stickyF(E, rho, kappa0)
%Evaluates sticky function and derivative for root finding. 
f=exp(E)/(sqrt(2*rho^2*E))-kappa0;
fprime=exp(E)/(2*rho^2*E)*(sqrt(2*rho^2*E)-rho^2/(sqrt(2*rho^2*E)));
end

