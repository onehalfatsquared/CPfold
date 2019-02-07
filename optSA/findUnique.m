function first=findUnique(sortedDistance,tol)
    %Takes a list of sorted particle distances and returns vector giving
    %the index where each unique distance first appears. Does so by
    %computing gaps between consecutive distances
    ele=length(sortedDistance); first(1)=1; %First always appears at 1
    for i=1:ele-1
        gap=sortedDistance(i+1)-sortedDistance(i);
        if gap>tol
            first=[first,i+1];
        end
    end
end