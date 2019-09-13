function u = morseEval(particles,rho,E,N,P)
    %evaluate morse potential for cluster
    u = 0; %initialize potential
    rep = 0; %repulsion coeff
    for i = 1:N
        for j = i+1:N
            %if particles are interacting, compute sum contribution
            if j~=i 
                %Morse interaction
                r = euDist(particles,i,j);
                if P(min(i,j),max(i,j))==1
                    u = u + morseD(r,rho,E(min(i,j),max(i,j))); 
                else
                    %steric repulsion
                    if r<1
                        u = u - rep/r;
                    end
                end
            end
        end
    end
end
