function l = fit_line_e(p1, p2, e)

p1 = p1-e;
p2 = p2-e;

A = [p1'; p2'];

% min norm(Ax) s.t. norm(x)=1
[~,~,V] = svd(A);

% normal
n = V(:,end);

% intercept
c = -n(1)*e(1)-n(2)*e(2);
l = [n; c];

end