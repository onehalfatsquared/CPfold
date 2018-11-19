function g = morseGrad(particles,rho,E,N,P)
    %evaluate gradient of morse potential for each particle
    g = zeros(3*N,1); %initialize storage for gradient
    rep = 1.0; %repulsion coeff
    for i = 1:N
        S = 0; %initialize sum
        for j = 1:N
            %if particles are interacting, compute sum contribution
            if j~=i 
                r = euDist(particles,i,j);
                if P(min(i,j),max(i,j))==1
                    S = S + morseD(r,rho,E)/r*(particles(i,:)-particles(j,:))';
                else
                    if r<1
                        S = S+rep*(1-1/r)*1/r^2;
                    end
                end
            end
        end
        g(3*(i-1)+1:3*i) = S;
    end
end