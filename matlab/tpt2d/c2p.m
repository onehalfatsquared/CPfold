function particles=c2p(cluster)
    %Seperate cluster vector into particle array
    n=length(cluster)/2;
    particles = reshape(cluster,2,n)';
end
