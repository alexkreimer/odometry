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

N = 1;
H_inf = nan(3,3,N);

R = estimation.H_inf_nonlin(K, x1, x2, d(inliers), b, thr1, thr2);
% initial computation
H_inf(:,:,1) = K*R/K;

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
T = nan(4,4,N);
for i=1:size(H_inf,3)
    H = H_inf(:,:,i);
    R = K\H*K;

    % warp the points. Note, we warp 2->1 with inv(H)
    x1ot = util.h2e(H*x1o);

    % find the epipole
    a = util.h2e(x2o);
    b = x1ot;

    lines = compute_initial_lines(a, b);
    e0 = util.h2e(cross(lines(:,1),lines(:,2)));

    opt = optimset(optimset('lsqnonlin') , 'Algorithm','levenberg-marquardt', 'Diagnostics','off', 'Display','off');
    fun = @(e) objective(a, b, e);
    [e,~,~,~] = lsqnonlin(fun, e0, [], [], opt);
    %if norm(e-e_gt)/norm(e_gt) > 1
    %    warn('epipole is way off');
    %end
    
    t = K\util.e2h(e);
    t = t/norm(t);
    T(:,:,i) = [R t; 0 0 0 1];
end

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

function val = objective(a, b, e)

% a 2xN points s.t. a(i) is a point on line i
% b 2xN points s.t. b(i) is a point on line i
% e epipole

N   = size(a,2);
val = zeros(2*N,1);

for i=1:N
    p1 = a(:,i);
    p2 = b(:,i);

    l = estimation.fit_line_e(p1, p2, e);
    
    val(2*i-1) = l'*[p1;1];
    val(2*i-0) = l'*[p2;1];
end
end

function lines = compute_initial_lines(a, b)
N     = size(a,2);
lines = nan(3,N);

for i = 1:N
    p1 = [a(:,i); 1];
    p2 = [b(:,i); 1];
    lines(:,i) = cross(p1, p2);
end
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

