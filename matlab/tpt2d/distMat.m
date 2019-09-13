function [D] = distMat(clust)
    %outputs distance matrix of a cluster
    particles = c2p(clust); L=length(particles); D = zeros(L);
    for i=1:L
        for j=i+1:L
            D(i,j)=euDist(particles,i,j);
        end
    end
end

