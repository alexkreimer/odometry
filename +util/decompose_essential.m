function T = decompose_essential(E, K, x)
% Decompose essential matrix E and resolve the twisted pair ambiguity by
% veryfing the depth of the reconstructed points
%
% K is a 3 by 4 intrinsics matrix
% x is 4 by N matrix. x(1:2,j)/x(3:4,j) is a left/right image projection
% for point j


[U,S,V] = svd(E);

% H&Z solution

% rotz(pi/2)
W = [0,-1,0;+1,0,0;0,0,1];
Z = [0,+1,0;-1,0,0;0,0,0];
t = util.vex(U*Z*U');
t = t/norm(t);

R1 = U*W *V';
if det(R1)<0
    R1 = -R1;
end

R2 = U*W'*V';
if det(R2)<0
    R2 = -R2;
end

% possible solutions
s(:,:,1) = [R1 +t];
s(:,:,2) = [R2 +t];
s(:,:,3) = [R1 -t];
s(:,:,4) = [R2 -t];

num_pts = size(x,2);
max_inliers = 0;
% reconstruct x
for i=1:size(s,3)
    P1 = K*[eye(3) zeros(3,1)];
    P2 = K*s(:,:,i);
    detM1 = det(P1(1:3,1:3));
    detM2 = det(P2(1:3,1:3));
    
    inliers = false(num_pts,1);
    for j=1:num_pts
        p = [x(1:2,j),x(3:4,j)];
    
        % Linear triangulation;
        % http://www.robots.ox.ac.uk/~vgg/hzbook/code/vgg_multiview/vgg_X_from_xP_lin.m
        X = vgg.vgg_X_from_xP_lin(p,{P1,P2});

        % enforce chierality constraint
        x1 = P1*X;
        x2 = P2*X;
        if x1(3)*X(4)*detM1>0 && x2(3)*X(4)*detM2>0
            inliers(j) = true;
        end
    end
    
    % majority vote.  their might be some outliers
    n = sum(inliers);
%     fprintf('total points %d, inliers %d\n', length(x), n);
    if n>max_inliers
        T = s(:,:,i);
        max_inliers = n;
    end
end
end
