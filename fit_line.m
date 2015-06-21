function [c,n] = fit_line(pts)

% pts is 2 \times n, where n is the number of points
% A (c n)' ~ 0 subject to norm(n,2)=1
% c is the intercept and n is the normal

x = pts(1,:)';
y = pts(2,:)';
A = [ones(3,1),x,y];
[c,n] = clsq(A,2);