function [E] = createE(types,K)
    %create matrix of well-depths specified by types of particles
    %interacting. Supports 2 types. A->A, A<->B, or B->B. 
    N = length(types); %Number of particles
    E = zeros(N);      %Initialize pairwise interaction matrix
    for i = 1:N
        for j = i+1:N
            if types(i) == 1 && types(j) == 1 %Both type A's
                E(i,j) = K(1);
            elseif types(i) == 2 && types(j) == 2 %Both types B's
                E(i,j) = K(3); 
            else                               %Both different
                E(i,j) = K(2);
            end
        end
    end
end