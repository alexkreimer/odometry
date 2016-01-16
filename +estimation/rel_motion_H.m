function [T,inliers] = rel_motion_H(K, x1, x2, d, b, thr1, thr2, e_gt)

% x1 interest points in image 1
% x2 interest points in image 2
% d  depth of the points
% b  baseline
% size(x1)=size(x2)=[2,N]

if nargin<7
    thr2 = .01;
end
if nargin<6
    thr1 = 250;
end

assert(all(size(x1)==size(x2)));

[rows, ~] = size(x1);

if rows==2
    x1 = util.e2h(x1);
    x2 = util.e2h(x2);
end

inliers = d>thr1*b;
if sum(inliers<10)
    inliers = d>.5*thr1*b;
    if sum(inliers<10)
        inliers = d>.5*thr1*b;
    end
end
x1o = x1(:,~inliers);
x2o = x2(:,~inliers);
x1  = x1(:, inliers);
x2  = x2(:, inliers);

R = estimation.H_inf_nonlin(K, x1, x2, d(inliers), b, thr1, thr2);
H = K*R/K;

% objective1
% H   = H(:,:,1)/H(3,3,1);
% h0  = H(1:8);
% fun = @(h) objective1(x1(:,inliers_H),x2(:,inliers_H),F,h);
% [h, resnorm, residual, exitflat] = lsqnonlin(fun, h0);
% H = ones(3);
% H(1:8) = h;
% H_inf(:,:,2) = H;
% H_gt = K*R_gt/K;

% figure;
% hold on;
% title('H error vs. point depth on log-scale');
% x1t_gt = util.h2e(H_gt*x1_orig);
% err_gt = sqrt(.5*sum((util.h2e(x2_orig)-x1t_gt).*(util.h2e(x2_orig)-x1t_gt)));
% semilogy(d,err_gt,'o','DisplayName','gt');
% names = {'H','H-objective1'};
% for j=1:size(H_inf,3)
%     x1t = util.h2e(H_inf(:,:,j)*x1_orig);
%     err = sqrt(.5*sum((util.h2e(x2_orig)-x1t).*(util.h2e(x2_orig)-x1t)));
%     semilogy(d,err,'o','DisplayName',names{j});
% end
% legend show;


t = estimation.trans_geom(K, H, x1o, x2o);
e = util.h2e(K*t);

if nargin >7 && norm(e-e_gt)/norm(e_gt) > 1
    warn('epipole is way off');
end

T = [R t; 0 0 0 1];

% figure; hold on;
% for i = 1:N
%      x_min = min([a(1,i) b(1,i) e(1)]);
%      x_max = max([a(1,i) b(1,i) e(1)]);
%      x = x_min:x_max;
%      l = estimation.fit_line_e(a(:,i),b(:,i),e);
%      y = util.get_line_points(l(1:2),e,x);
%      plot(x,y);
%      plot(a(1,i), a(2,i),'og');
%      plot(b(1,i), b(2,i),'ob');    
% end
% 
% plot(e(1),e(2),'*r', 'DisplayName', 'estimated epipole');
% plot(e_gt(1), e_gt(2), '*g', 'DisplayName', 'gt epipole');
% legend show;

% OUT_DIR = '/home/kreimer/prj/odometry/debug';
% save_dbg(fullfile(OUT_DIR, tempname));
% close;
end



function val = objective1(x1,x2,F,h)
% symmetric transfer error
w1 = .8;
w2 = 1;
val = nan(1,2*length(x1));

H = ones(3);
H(1:8) = h;

for i = 1:length(x2)
    p2 = x2(:,i);
    p1 = x1(:,i);
    
    delta1 = util.h2e(H*p1)-util.h2e(p2);
    delta2 = util.h2e(H\p2)-util.h2e(p1);
    
    l1 = F*p1;
    n1 = l1(1:2)/norm(l1(1:2));
    v1 = [0 -1;1 0]*n1;
    
    l2 = F'*p2;
    n2 = l2(1:2)/norm(l2(1:2));
    v2 = [0 -1;1 0]*n2;
    
    val(2*i-1) = w1*n1'*delta1 + w2*v1'*delta1;
    val(2*i)   = w1*n2'*delta2 + w2*v2'*delta2;
end
end

