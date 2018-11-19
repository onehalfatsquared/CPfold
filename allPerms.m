function AP=allPerms(f,r)
    %Construct all possible permutations of particles consistent with their
    %distances from the origin.
    
    %Create the grouping from f and r
    grouping={1};     %Initialize grouping
    n=length(r);      %Number of unique distances
    if r(1)==1
        for i=2:n
            grouping{i}=f(i):f(i)+r(i)-1;
        end
    else
        grouping{2}=2:f(1)+r(1)-1;
        for i=3:n+1
            grouping{i}=f(i-1):f(i-1)+r(i-1)-1;
        end
    end
    
    %Create the array of permutations consistent with grouping
    AP=[];     %Initialize array
    perm=cellfun(@(x){perms(x)},grouping); %for each group generate all permutations
    cart=cell(1,numel(grouping));  %Initialize cartesian product storage
    perm_size=cellfun(@(x){1:size(x,1)},perm); %[1:size_of_permutation] for ndgrid
    [cart{:}]=ndgrid(perm_size{:});  %Cartesian product of indexes of permutations
    for k=1:numel(cart) %index permutations with the cartesian product of their indexes
        AP=horzcat(AP,perm{k}(cart{k},:));
    end
end