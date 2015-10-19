function val = sampson(q,t,x1,x2)

% q - unit quaternion that represents rotation
% t - translation
% E = [t]xR
% x1,x2 inlier coordinates, s.t. x2'Ex1=0

N = size(x1,2);

R = quaternion(q).RotationMatrix;

E = skew(t)*R;
val = nan(N,1);
for i=1:N
    % current point
    p1 = x1(:,i);
    p2 = x2(:,i);
    
    p2tE = p2'*E;
    Ep1  = E*p1;
    
    val(i) = p2'*skew(t)*R*p1/sqrt(p2tE(1)^2 + (p2tE(1))^2 + (Ep1(1))^2 + (Ep1(2))^2);
end
end


