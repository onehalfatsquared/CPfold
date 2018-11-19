function U=findRot(P,Q)
    %Find optimal rotation to bring set of particles P to set of particles
    %Q. Kabsch algorithm
    A=P'*Q;             %Covariance matrix
    [V,S,W]=svd(A);     %SVD of A 
    %d=det(W*V');       %Determines handedness of coordinate system. Need right
    d=1;                %d=1 gives possibility of reflection+rotation
    U=W*diag([1,1,d])*V'; %Optimal rotation matrix
end