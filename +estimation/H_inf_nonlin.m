function [H,R_out,inliers_out] = H_inf_nonlin(K,F,x1,x2,R0)

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

% if nargin > 4
%     H0 = K*R0/K;
%     H0x1 = util.h2e(H0*x1);
% end
% 
% % The 3x3 homography such that x2 = H*x1.
% H = vgg.vgg_H_from_x_lin(x1,x2);
% [~,d] = util.homogdist2d(H,x1,x2,0);
% fprintf('homography reprojection error is RMS is %g\n',sqrt(d*d'/length(d)));
% 
% Hx1 = util.h2e(H*x1);
% 
% figure;
% hold on;
% 
% Fx1  = F*x1;
% Fx1  = Fx1./repmat(util.colnorm(Fx1(1:2,:)),[3 1]);
% l1   = Fx1(1:2,:)'*[0 1;-1 0];
% l1x2 = util.h2e(x2) - l1';
% 
% plot([x2(1,:); Hx1(1,:)], [x2(2,:); Hx1(2,:)], '-og');
% plot([x2(1,:); l1x2(1,:)],[x2(2,:); l1x2(2,:)],'-oc');
% 
% if nargin>4
%     plot([x2(1,:); H0x1(1,:)],[x2(2,:); H0x1(2,:)],'-ob');
% end
% 
%     R = K\H*K;
%     [U,S,V] = svd(R)
%     R = U*eye(3)*V';
%     H = K*R/K;
%
%     Hx1 = util.h2e(H*x1);
%     hold on;
%     plot([x2(1,:); Hx1(1,:)],[x2(2,:); Hx1(2,:)],'-or', 'DisplayName', 'Normalized homography');
%
%     [~,d] = util.homogdist2d(H,x1,x2,0);
%     fprintf('homography forced to KR/K reprojection error is RMS is %g\n',sqrt(d*d'/length(d)));
% r = vrrotmat2vec(R);
% vt0 = r(1:3)*r(4);
num_pts = length(x1);
s = 4;
max_size = 0;
thr2 = 2;
for i = 1:100
    sample = randsample(num_pts,s);
    vt0 = zeros(1,3);
    fun = @(vt) objective1(K,F,x1(:,sample),x2(:,sample),vt);
    [vt,resnorm,result,exitflag] = lsqnonlin(fun,vt0);
    R = vrrotvec2mat(vtheta2r(vt));
    H = K*R/K;
    [~,d] = util.homogdist2d(H,x1(:,sample),x2(:,sample),0);
    fprintf('homography symm reprojection error RMS is %g\n',sqrt(d*d'/length(d)));
    inliers = util.homogdist2d(H,x1,x2,thr2);
    support_size = sum(inliers);
    fprintf('support size: %d\n',support_size);
    if support_size>max_size
        max_size = support_size;
        inliers_best = inliers;
        vt_best = vt;
    end
end

fun = @(vt) objective1(K,F,x1(:,inliers_best),x2(:,inliers_best),vt);
[vt,resnorm,result,exitflag] = lsqnonlin(fun,vt_best);
R = vrrotvec2mat(vtheta2r(vt));
H = K*R/K;
[~,d] = util.homogdist2d(H,x1(:,sample),x2(:,sample),0);
fprintf('homography symm reprojection error RMS is %g\n',sqrt(d*d'/length(d)));

% Hx1 = util.h2e(H*x1);
% plot([x2(1,:); Hx1(1,:)],[x2(2,:); Hx1(2,:)],'-om', 'DisplayName', 'Refined homography');
% legend show;
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

function val = objective1(K,F,x1,x2,vtheta)
r = vtheta2r(vtheta);
R = vrrotvec2mat(r);
H = K*R/K;
e1 = util.h2e(x2)-util.h2e(H*x1);
l1 = F*x1;
l1 = l1./repmat(util.colnorm(l1(1:2,:)),[3 1]);
e1_ortho = diag(l1(1:2,:)'*e1);

e2 = util.h2e(x1)-util.h2e(H\x2);
l2 = F'*x2;
l2 = l2./repmat(util.colnorm(l2(1:2,:)),[3 1]);
e2_ortho = diag(l2(1:2,:)'*e2);

val = e1_ortho.*e1_ortho + e2_ortho.*e2_ortho;
end