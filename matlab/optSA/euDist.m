function d=euDist(particles,i,j)
    %Compute the euclidean distance between particles i and j
    d=particles(j,:)-particles(i,:);
    %d=DNorm2(d,2);
    d=norm(d);
    %d=sqrt(d*d');
end