function [e, R] = epipole_from_H(K, x1, x2, d, b)

% x1 interest points in image 1
% x2 interest points in image 2
% d  depth of the points
% b  baseline
% size(x1)=size(x2)=[2,N]

assert(all(size(x1)==size(x2)));

[rows, ~] = size(x1);

if rows==2
    x1 = e2h(x1);
    x2 = e2h(x2);
end

inliers = d>30*b;
x1o = x1(:,~inliers);
x2o = x2(:,~inliers);
x1  = x1(:, inliers);
x2  = x2(:, inliers);

% The 3x3 homography such that x2 = H*x1.
[H, ~] = ransacfithomography_vgg(x1,x2,0.01);
R = K\H*K;

% warp the points
x2ot = h2e(H*x2o);
N    = size(x2ot,2);

% find the epipole
a = h2e(x1o);
b = x2ot;

lines = compute_initial_lines(a, b);
e0 = h2e(cross(lines(:,1),lines(:,2)));

opt = optimset( optimset('lsqnonlin') , 'Algorithm','levenberg-marquardt', 'Diagnostics','off', 'Display','off');
fun = @(e) mle_objective1(a, b, e);
[e,~,~,~] = lsqnonlin(fun, e0, [], [], opt);

% figure; hold on;
% for i = 1:N
% 
%     x_min = min([a(1,i) b(1,i) e(1)]);
%     x_max = max([a(1,i) b(1,i) e(1)]);
%     x = x_min:x_max;
%     l = fit_line_e(a(:,i),b(:,i),e);
%     y = get_line_points(l(1:2),e,x);
%     plot(x,y);
%     plot(a(1,i), a(2,i),'og');
%     plot(b(1,i), b(2,i),'ob');    
% end
% 
% plot(e(1),e(2),'*r');

% OUT_DIR = '/home/kreimer/prj/odometry/debug';
% save_dbg(fullfile(OUT_DIR, tempname));
% close;
end

function val = mle_objective(a, b, n, v)

% argmin 1/N( sum_{1}^{N} d^2(l(i),a(i)) + d^2(l(i),b(i)) 
%   s.t. 1/N(sum_{1}^{N} d^2(l(i),v) = 0

% a 2xN points s.t. a(i) is a point on line i
% b 2xN points s.t. b(i) is a point on line i
% v is the itersection point
% n are the normals

N   = size(a,2);
val = zeros(2*N+N,1);

for i=1:N
    % current normal vector
    n_cur = n((2*i-1):2*i);
    % current line equation
    l_cur = [n_cur' -n_cur'*v]';
    
    val(2*i-1) = line_ortho_dist(l_cur, a(:,i));
    val(2*i-0) = line_ortho_dist(l_cur, b(:,i));
    
    val(2*N+i) = line_ortho_dist(l_cur,v);
end
end

function val = mle_objective1(a, b, e)

% a 2xN points s.t. a(i) is a point on line i
% b 2xN points s.t. b(i) is a point on line i
% e epipole

N   = size(a,2);
val = zeros(2*N,1);

for i=1:N
    p1 = a(:,i);
    p2 = b(:,i);

    l = fit_line_e(p1, p2, e);
    
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


