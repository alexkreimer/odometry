function e = epipole_from_H(x1, x2)

% x1 interest points in image 1
% x2 interest points in image 2
%
% size(x1)=size(x2)=[2,N]

assert(all(size(x1)==size(x2)));

[rows, cols] = size(x1);

if rows==2
    x1 = e2h(x1);
    x2 = e2h(x2);
end

% The 3x3 homography such that x2 = H*x1.
[H, inliers] = ransacfithomography_vgg(x1,x2,0.0001);

% figure; hold on;
% plot([x1(1,:); x2(1,:)], [x1(2,:); x2(2,:)]);
x1_inliers = h2e(x1(:, inliers));
x2_inliers = h2e(x2(:, inliers));

% warped points
x2t = h2e(H*x1);
x2t_inliers = x2t(:,inliers);

% h1 = plot(x1_inliers(1,:), x1_inliers(2,:), 'og');
% h2 = plot(x2t(1,:), x2t(2,:), 'ob');
% h3 = plot(x2t_inliers(1,:), x2t_inliers(2,:), 'or');
% legend([h1,h2,h3],'left inliers', 'all left points warped', 'inliers warped');

outliers = true([1,size(x1,2)]);
outliers(inliers) = false;

N = size(x1,2);
for i=1:N
    if norm(x2t(:,i)-h2e(x1(:,i)))<2
        outliers(i) = true;
    end
end

x2t_out = x2t(:, outliers);
x1_out = h2e(x1(:, outliers));

% figure; hold on; title('remaining interest points after H was found');
% plot(x2t_out(1,:), x2t_out(2,:), 'og');
% plot(x1_out(1,:), x1_out(2,:), 'or');
% legend('warped right', 'left');

N = size(x1_out,2);

lines = compute_initial_lines(x1_out, x2t_out);
v0 = h2e(cross(lines(:,1),lines(:,2)));
p0 = [lines(:); v0];

opt = optimset( optimset('lsqnonlin') , 'Algorithm','levenberg-marquardt', 'Diagnostics','off', 'Display','off');
fun = @(p) mle_objective(x1_out, x2t_out, p(1:end-2), p(end-1:end));
[p,resnorm,residual,exitflag] = lsqnonlin(fun, p0, [], [], opt);

e = p(end-1:end);

for i = 1:N
    x_min = min([x1_out(1,i) x2t_out(1,i) e(1)]);
    x_max = max([x1_out(1,i) x2t_out(1,i) e(1)]);
    
    x = x_min:x_max;
    y = get_line_points(p((2*i-1):2*i),e,x);
    plot(x,y);
end

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

function d = line_ortho_dist(l,x)
% distance from point x to line l

n = (l(1)*l(1)+l(2)*l(2));
d = (l(1)*x(1)+l(2)*x(2)+l(3))/n;
   
d = abs(d);
end

function lines = compute_initial_lines(a, b)

N = size(a,2);
lines = nan(3,N);

for i = 1:N
    p1 = [a(:,i);1];
    p2 = [b(:,i);1];
    
    lines(:,i) = cross(p1,p2);
end

end

function y = get_line_points(n, v, x)

a = -n(1)/n(2);
b = n'*v/n(2);
y = a*x + b;


end
