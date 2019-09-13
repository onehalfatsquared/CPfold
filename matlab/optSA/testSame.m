function [b,Pout,Qout,rmsd,L1diff,minDist]=testSame(c1,c2,arg3)
    %Test if clusters 1 and 2 are the same up to translation, rotation, and
    %permutation. Optionally, determine distance between clusters if arg3 given.  
    
    %Argument management. arg3 is optional. If provided and not 0, set to 1.
    if nargin<3 || arg3==0
        arg3=0;
    else
        arg3=1;
    end
    
    %No arg3 -> set parameters just to test if clusters are the same. 
    if arg3==0
        b=0; Pout=NaN; Qout=NaN; rmsd=Inf; L1diff=Inf; minDist=Inf; %Initialize with output for not the same cluster
        potTol=1e-7; distTol=1e-7; sameTol=1e-5; gapTol=0.1; %Tolerances
    end
    
    %arg3 given -> set parameters to test if same, but record minimum
    %distance between them if they are different. 
    if arg3==1
        b=0; Pout=NaN; Qout=NaN; rmsd=Inf; L1diff=Inf; minDist=Inf; %Initialize with output for finding min
        potTol=Inf; distTol=Inf; sameTol=1e-6; gapTol=1; %Tolerances
    end
    
    %Convert clusters to particle arrays. 
    p1=c2p(c1); p2=c2p(c2);
    
    %First compute energy to see if clusters are potentially the same
    rho=30; E=30;
    E1=MP(p1,rho,E); E2=MP(p2,rho,E);
    if abs(E1-E2)>potTol
        return
    end
    
    %Subtract the "Center of Mass" from each set of particles
    p1=p1-ones(size(p1))*diag(sum(p1)/length(p1));
    p2=p2-ones(size(p2))*diag(sum(p2)/length(p2));
    
    %Compute distances to origin of each particle
    d1=sqrt(sum(p1.^2,2));
    d2=sqrt(sum(p2.^2,2));
    
    %Sort these lists of distances. Compare within tol. If same, permute
    %particle lists
    [sd1,I1]=sort(d1); [sd2,I2]=sort(d2);
    if max(abs(sd1-sd2))>distTol
        return
    end
    sp1=p1(I1,:); sp2=p2(I2,:);
    
    %Count the unique distances, determine how many times each occurs
    first=findUnique(sd1,gapTol); %Where each distance first appears
    for i=1:length(first)-1
        times(i)=first(i+1)-first(i);
    end
    times(length(first))=length(d1)+1-first(length(first)); %Num times appears
    
    %Choose first index of least appearing distance. How many repeats?
    indexL=first(find(times==min(times),1)); %Least appearing index
    repeatL=times(ismember(first,indexL)==1); %Number of times appears
    
    %Pick first index. Rotate this particle to all particles in second
    %cluster with the same distance. Check if same cluster, if not permute
    u=sp1(indexL,:);          %Particle being rotated
    AP=allPerms(first,times); %All consistent perms of particles
    for test=1:repeatL
        P=sp1; Q=sp2;    %Collection of particles being rotated/ rotated to
        Q([1,test],:)=Q([test,1],:); %Swap particle being rotated to every iter.
        v=Q(indexL,:);   %Candidate particle to align with
        R=vecRot(u,v);   %Find rotation to take u to v
        
        for i=1:length(P)
            P(i,:)=P(i,:)*R'; %Apply rotation to each particle
        end
        
        for i=1:size(AP,1)
            Ptest=P(AP(i,:),:);     %Permutation being tested
            U=findRot(Ptest,Q);     %Optimal rotation of Ptest to sp2
            for j=1:length(Ptest)
                Ptest(j,:)=Ptest(j,:)*U'; %Apply rotation to each particle
            end
            rmsd=sqrt(1/length(Ptest)*sum(sum((Ptest-Q).^2)));    %L2 distance
            if rmsd<minDist
                minDist=rmsd; %Min distance among perms
                Pout=Ptest; Qout=Q; %Closest configurations
            end
            L1diff=sum(sum(abs(Ptest-Q))); %L1 distance
            if rmsd<sameTol
                b=1;       %Same cluster
                minDist=rmsd;
                return
            end
        end  
    end
end