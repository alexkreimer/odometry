function val = sampson_err(R,t,x1,x2)

% R - rotation matrix
% t - translation
% E = [t]xR
% x1,x2 inlier coordinates, s.t. x2'Ex1=0

N = size(x1,2);

E = skew(t)*R;
val = nan(N,1);
for i=1:N
    % current point
    p1 = x1(:,i);
    p2 = x2(:,i);
    
    if length(p1) == 2
        p1 = [p1;1];
    end
    
    if length(p2) == 2
        p2 = [p2;1];
    end
    
    
    p2tE = p2'*E;
    Ep1  = E*p1;
    
    val(i) = p2'*skew(t)*R*p1/sqrt(p2tE(1)^2 + (p2tE(1))^2 + (Ep1(1))^2 + (Ep1(2))^2);
end
end