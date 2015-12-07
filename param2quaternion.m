function hZ = param2quaternion(v, h0)
% convert 3-vector v and the operating point h0 into a quaternion hZ

B = ortho_basis(h0);

% convert parameter vector to 4-vector and normalize it
v4  = B*v;

if norm(v4) == 0
    hZ = h0;
else
    v4N = v4/norm(v4);
    
    % resulting quaternion
    theta = norm(v4);
    hZ = sin(theta)*v4N + cos(theta)*h0;
end

end

function B = ortho_basis(h)

if h(1) ~= 0
    b0 = [1/h(1); 0; 0; 0];
    b1 = [-h(2)/h(1); 1; 0; 0];
    b2 = [-h(3)/h(1); 0; 1; 0];
    b3 = [-h(4)/h(1); 0; 0; 1];
    
    % non orthogonal basis
    B = [b1 b2 b3];
    
    [U,~,~] = svd(B);
    
    % orthogonal basis that spans the same subspace
    B = U(:,1:3);
end

end
