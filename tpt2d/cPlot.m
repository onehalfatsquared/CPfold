function [] = cPlot(x,types,P)
    %Plot cluster of particles with edges
    
    %create color scheme for particle types
    N = length(types); color = zeros(N,3);
    for i=1:length(types)
        if types(i)==1
            color(i,:)=[1,0,0];
        elseif types(i)==2
            color(i,:)=[0,1,0];
        elseif types(i)==3
            color(i,:)=[0,0,1];
        end
    end
    
    %add 0 z components to easily plot
    temp = reshape(x, 2, N);
    temp = temp'; temp(:,3) = zeros(N,1);
    temp = temp'; temp = reshape(temp, 1, 3*N);
    x = temp';
    
            
    %matrix to kill bonds for non-interacting particles
    Z = P+P';
    
    %options for plotting cluster
    opts = struct('srad',0.1,'lcolr',0.7*[1,1,1],'lrad',0.03,...
    'scolr',color,'salph',0.7,'ifgrid',0,'ax',[-N,N,-N/2,N/2,-N/3,N/3]/2,...
    'iftight',0,'az',20,'el',18,'zth',-0.5,'special',Z);

    %plot the cluster
    plotcluster2b(x,opts); 
end

