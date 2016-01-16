function R = H_inf(K, x1, x2, d, b, thr1, thr2)

% x1 interest points in image 1
% x2 interest points in image 2
% d  depth of the points
% b  baseline
% size(x1)=size(x2)=[2,N]

assert(all(size(x1)==size(x2)));

[rows, ~] = size(x1);

if rows==2
    x1 = util.e2h(x1);
    x2 = util.e2h(x2);
end

inliers = d>thr1*b;

x1  = x1(:, inliers);
x2  = x2(:, inliers);

% The 3x3 homography such that x2 = H*x1.
[H, inliers_H] = estimation.ransacfithomography_vgg(x1, x2, thr2);

R = K\H*K;
[u,s,v] = svd(R);
s
R = u*eye(3)*v';
end