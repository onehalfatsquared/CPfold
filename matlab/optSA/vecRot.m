function R=vecRot(u,v)
    %Find rotation matrix in R^3 that takes u to v. 
    u=u/norm(u); v=v/norm(v);
    ssc = @(vec) [0 -vec(3) vec(2); vec(3) 0 -vec(1); -vec(2) vec(1) 0];
    R=eye(3)+ssc(cross(u,v))+ssc(cross(u,v))^2*(1-dot(u,v))/(norm(cross(u,v))^2);
end