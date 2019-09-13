function Ur = morseD(r,rho,E)
    %Evaluate morse potential derivative
    Y=exp(-rho*(r-1));
    Ur=-2*rho*E*(Y^2-Y);
end