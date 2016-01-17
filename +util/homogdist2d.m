function [inliers, d] = homogdist2d(H,x1,x2,t)
% x2 = Hx1

% Calculate, in both directions, the transfered points
Hx1    = H*x1;
invHx2 = H\x2;

x1     = util.h2e(x1);
x2     = util.h2e(x2);
Hx1    = util.h2e(Hx1);
invHx2 = util.h2e(invHx2);

d2 = sum((x1-invHx2).^2)  + sum((x2-Hx1).^2);
inliers = find(abs(d2) < 4*t*t);

if nargout==2
    d = d2;
end

end