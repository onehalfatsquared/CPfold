function [E] = createE(types,Ek)
    %create matrix of well-depths specified by types of particles
    %interacting. Supports 2 types. A->A, A<->B, or B->B. 
    N = length(types); %Number of particles
    E = zeros(N);      %Initialize pairwise interaction matrix
    for i = 1:N
        for j = i+1:N
            if types(i) == 1 && types(j) == 1 %Both type A's
                E(i,j) = Ek(1);
            elseif types(i) == 2 && types(j) == 2 %Both types B's
                E(i,j) = Ek(3); 
            else                               %Both different
                E(i,j) = Ek(2);
            end
        end
    end
end