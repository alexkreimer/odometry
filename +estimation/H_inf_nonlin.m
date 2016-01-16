function R = H_inf_nonlin(K, x1, x2, d, b, thr1, thr2)

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
[H, inliers_H] = estimation.ransacfithomography_vgg(x1, x2, thr2, 1);

R = K\H*K;
[U,S,V] = svd(R);
R = U*eye(3)*V';
r = vrrotmat2vec(R);
vt0 = r(1:3)*r(4);
fun = @(vt) objective1(K,x1,x2,vt);
[vt,resnorm,result,exitflag] = lsqnonlin(fun,vt0);
R = vrrotvec2mat(vtheta2r(vt));
end

function r = vtheta2r(vtheta)
    angle = norm(vtheta);
    if angle == 0
        axis = [0 1 0];
    else
        axis  = vtheta(1:3)/angle;
    end
    fprintf('axis: %f %f %f, angle %f\n', axis, angle);
    r = [axis angle];
end

function val = objective1(K,x1,x2,vtheta)
r = vtheta2r(vtheta);
R = vrrotvec2mat(r);
H = K*R/K;
val = util.h2e(x2)-util.h2e(H*x1) + util.h2e(x1)-util.h2e(H\x2);
val = sum(val.*val);
end