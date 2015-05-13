function T = decompose_essential(E,K,x)
% Decompose essential matrix E and resolve the twisted pair ambiguity by
% veryfing the depth of the reconstructed points
%
% K is a 3 by 4 intrinsics matrix
% x is 4 by N matrix. x(1:2,j)/x(3:4,j) is a left/right image projection
% for point j


[U,~,V] = svd(E);

% H&Z solution

% rotz(pi/2)
W = [0,-1,0;+1,0,0;0,0,1];
Z = [0,+1,0;-1,0,0;0,0,0];
t = vex(U*Z*U');
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
    % first camera
    P1 = K*[eye(3) zeros(3,1)];
    
    % 2nd camera. s transforms from the 1st camera into 2nd.

    P2 = K*s(:,:,i);
    
    P = nan(4,num_pts);
    inliers = false(num_pts,1);
    for j=1:num_pts
        p = [x(1:2,j),x(3:4,j)];
    
        % Linear triangulation;
        % http://www.robots.ox.ac.uk/~vgg/hzbook/code/vgg_multiview/vgg_X_from_xP_lin.m
        P(:,j) = vgg_X_from_xP_lin(p,{P1,P2});

        % enforce chierality constraint
        x2 = P2*P(:,j);
        if P(3,j)*P(4,j)>0 && x2(3)*P(4,j)>0
            inliers(j) = true;
        end
    end
    
    % majority vote.  their might be some outliers
    n = sum(inliers);
    if n>max_inliers
        T = s(:,:,i);
        max_inliers = n;
    end
end
end
