function particles=c2p(cluster)
    %Seperate cluster vector into particle array
    n=length(cluster)/3;
    particles=zeros(n,3);             %Store particle locations of cluster
    for i=1:n
        particles(i,:)=cluster(3*i-2:3*i);
    end
end