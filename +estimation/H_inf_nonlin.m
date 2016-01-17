function [H,R_out,inliers_out] = H_inf_nonlin(K, x1, x2, R0)

% x1 distant interest points in image 1
% x2 distant interest points in image 2
% d  depth of the points
% b  baseline
% size(x1)=size(x2)=[2,N]

assert(all(size(x1)==size(x2)));

[rows, ~] = size(x1);

if rows==2
    x1 = util.e2h(x1);
    x2 = util.e2h(x2);
end

if nargin < 8
    % The 3x3 homography such that x2 = H*x1.
    H = vgg.vgg_H_from_x_lin(x1,x2);
    
    R = K\H*K;
    [U,~,V] = svd(R);
    R = U*eye(3)*V';
else
    R = R0;
end

r = vrrotmat2vec(R);
vt0 = r(1:3)*r(4);
fun = @(vt) objective1(K,x1,x2,vt);
[vt,resnorm,result,exitflag] = lsqnonlin(fun,vt0);
R = vrrotvec2mat(vtheta2r(vt));

if nargout>1
    R_out = R;
    if nargout > 2
        inliers_out = inliers;
    end
end
end

function r = vtheta2r(vtheta)
angle = norm(vtheta);
if angle == 0
    axis = [0 1 0];
else
    axis  = vtheta(1:3)/angle;
end
r = [axis angle];
end

function val = objective1(K,x1,x2,vtheta)
r = vtheta2r(vtheta);
R = vrrotvec2mat(r);
H = K*R/K;
val = util.h2e(x2)-util.h2e(H*x1) + util.h2e(x1)-util.h2e(H\x2);
val = sum(val.*val);
end