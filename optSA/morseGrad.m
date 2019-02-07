function g = morseGrad(particles,rho,E,N,P)
    %evaluate gradient of morse potential for each particle
    g = zeros(3*N,1); %initialize storage for gradient
    rep = 3; %repulsion coeff
    for i = 1:N
        S = 0; %initialize sum
        for j = 1:N
            %if particles are interacting, compute sum contribution
            if j~=i 
                %Morse interaction
                r = euDist(particles,i,j);
                if P(min(i,j),max(i,j))==1
                    S = S + morseD(r,rho,E)/r*(particles(i,:)-particles(j,:))';
                else
                    %steric repulsion
                    if r<1
                        S = S - rep*(particles(i,:)-particles(j,:))'/r^3; 
                    end
                end
            end
        end
        g(3*(i-1)+1:3*i) = S;
    end
end