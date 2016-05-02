function t = trans_geom(K, H, x1, x2)
% find epipole by fitting a pencil of epipolar lines throught the points
% this code assumes pure camera translation, s.t. the features 'slide'
% along the epipolar lines
%
% x1 - feature points in the 1st view
% x2 - feature points in the 2nd view
% H  - infinite homography s.t. x2 = H*x1
%
% the optimization is over the epipole location; it minimizes a vector of
% geometric distances of the points to the epipolar lines

% warp the points

if size(x1,1) == 2
    x1 = util.e2h(x1);
    x2 = util.e2h(x2);
end
N   = length(x1);
Hx1 = util.h2e(H*x1);
a   = util.h2e(x2);
b   = Hx1;
thr = 2;
max_support = 0;
k = 0; kk = 0;
while k<5
    sample = randsample(N,2);
    p1 = a(:,sample(1));
    p2 = b(:,sample(1));
    l1 = cross(util.e2h(p1),util.e2h(p2));
    l1 = l1/norm(l1(1:2));
    p3 = a(:,sample(2));
    p4 = b(:,sample(2));
    l2 = cross(util.e2h(p3),util.e2h(p4));
    l2 = l2/norm(l2(1:2));
    if kk < 1000
        angle = acos(l1(1:2)'*l2(1:2));
        if abs(pi/2-angle)> pi/4
            kk = kk + 1;
            continue;
        end
    end
    k = k+1;
%     figure;
%     plot([p1(1); p2(1)],[p1(2); p2(2)],'.b');
%     hold on;
%     plot([p3(1);p4(1)],[p3(2);p4(2)],'.r');
    
    % epipole hypothesis
    e = util.h2e(cross(l1, l2));
    
%     x_min = min([p1(1) p2(1) e(1)]);
%     x_max = max([p1(1) p2(1) e(1)]);
%     x = x_min:x_max;
%     y1 = util.get_line_points(l1(1:2),e,x);
%     plot(x,y1,'DisplayName', 'hypothesis line1');
%     
%     x_min = min([p3(1) p4(1) e(1)]);
%     x_max = max([p3(1) p4(1) e(1)]);
%     x = x_min:x_max;
%     y2 = util.get_line_points(l2(1:2),e,x);
%     plot(x,y2,'DisplayName', 'hypothesis line2');
%     plot(e(1),e(2),'*r', 'DisplayName', 'epipole hypothesis');

    inliers = false(1,N);
    for j = 1:N
%         x_min = min([a(1,j) b(1,j) e(1)]);
%         x_max = max([a(1,j) b(1,j) e(1)]);
%         x = x_min:x_max;
        
        [l,d1,d2] = estimation.fit_line_e(a(:,j),b(:,j),e);
%         y = util.get_line_points(l(1:2),e,x);
        
        inliers(j) = abs(d1)<thr && abs(d2)<thr;
%         if inliers(j)
%             plot(x,y,'g--');
%             plot(a(1,j), a(2,j),'og');
%             plot(b(1,j), b(2,j),'ob');
%         else
%             plot(x,y,'y--');
%             plot(a(1,j), a(2,j),'*g');
%             plot(b(1,j), b(2,j),'*b');
%         end
    end
    n = sum(inliers);
    if n>max_support
        best_inliers = inliers;
        best_e = e;
        max_support = n;
    end
end

opt = optimset(optimset('lsqnonlin') , 'Algorithm','levenberg-marquardt', 'Diagnostics','off', 'Display','off');
fun = @(e) objective(a(:,best_inliers), b(:,best_inliers), e);
[e,resnorm,residual,exitflag] = lsqnonlin(fun, best_e, [], [], opt);
t = K\[e;1];
t = t/norm(t);

end

function lines = compute_initial_lines(a, b)
% intersection of 2 first epipolar lines
N     = size(a,2);
lines = nan(3,N);

for i = 1:N
    p1 = [a(:,i); 1];
    p2 = [b(:,i); 1];
    lines(:,i) = cross(p1, p2);
end
end

function res = objective(a, b, e)
% a 2xN points s.t. a(i) is a point on line i
% b 2xN points s.t. b(i) is a point on line i
% e epipole
N   = size(a,2);
res = zeros(2*N,1);

for i=1:N
    p1 = a(:,i);
    p2 = b(:,i);
    
    l = estimation.fit_line_e(p1, p2, e);
    
    res(2*i-1) = l'*[p1;1];
    res(2*i-0) = l'*[p2;1];
end

if 0
    % reweight
    weights = exp(-3*abs(res)/max(abs(res)));
    res = res.*weights;
end
end

function [n,inliers] = compute_support(x1,x2,e,thr)

N = length(x1);
inliers = false(1,N);
for i=1:N
    p1 = x1(:,i)-e;
    p2 = x2(:,i)-e;
    
    A = [p1'; p2'];
    
    % min norm(Ax) s.t. norm(x)=1
    [~,~,V] = svd(A);
    
    % normal
    n = V(:,end);
    
    % intercept
    c = -n(1)*e(1)-n(2)*e(2);
    l = [n; c];
    l = l/norm(l(1:2));
    
    d1 = l'*util.e2h(x1(:,i));
    d2 = l'*util.e2h(x2(:,i));
    
    inliers(i) = d1 < thr && d2 < thr;
end
n = sum(inliers);
end