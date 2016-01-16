function t = trans_geom(K,H,x1,x2)
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
x1t = util.h2e(H*x1);

a = util.h2e(x2);
b = x1t;

lines = compute_initial_lines(a, b);
e0 = util.h2e(cross(lines(:,1),lines(:,2)));

opt = optimset(optimset('lsqnonlin') , 'Algorithm','levenberg-marquardt', 'Diagnostics','off', 'Display','off');
fun = @(e) objective(a, b, e);
[e,resnorm,residual,exitflag] = lsqnonlin(fun, e0, [], [], opt);

t = K\e;
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
