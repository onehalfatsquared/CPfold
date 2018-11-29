function particles=c2p(cluster)
    %Seperate cluster vector into particle array
    n=length(cluster)/3;
    particles = reshape(cluster,3,n)';
end