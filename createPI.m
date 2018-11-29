function [P] = createPI(types)
    %create pairwise interaction matrix for each particle
    %type 1 bonds to everything, type 2 and 3 do not bind eachother
    N = length(types); %Number of particles
    P = zeros(N);      %Initialize pairwise interaction matrix
    for i=1:N
        for j=i+1:N
            if types(i)*types(j) ~=6
                P(i,j) = 1;
            end
        end
    end   
end