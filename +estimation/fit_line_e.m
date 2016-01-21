function [l,d1,d2] = fit_line_e(p1, p2, e)

p1e = p1-e;
p2e = p2-e;

A = [p1e'; p2e'];

% min norm(Ax) s.t. norm(x)=1
[~,~,V] = svd(A);

% normal
n = V(:,end);

% intercept
c = -n(1)*e(1)-n(2)*e(2);
l = [n; c];

if nargout>1
    d1 = l'*util.e2h(p1);
    d2 = l'*util.e2h(p2);
end

end