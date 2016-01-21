function simple_H_inf()
close all;
dbstop if error;

K = [1 0 0; 0 1 0; 0 0 1];
% K = [718.8560,  0,      607.1928;
%        0,     718.8560, 185.2157;
%        0,       0,        1.0000];       % camera intrinsics
   
N = 100;
X = util.e2h(10+10*randn(3,N));

I = eye(3);
%R = rotx(randn*.2)*roty(randn*.2);
R = I;
t = [0 0 2]';

x1 = K*[I zeros(3,1)]*X;
x2 = K*[R' -R'*t]*X;


H = K*R'/K;

Hx1   = H*x1;
delta = x2-Hx1;

figure;
p1 = util.h2e(x1);
p2 = util.h2e(Hx1);
plot([p2(1,:); p1(1,:)],[p2(2,:); p1(2,:)]);

% sample mean
m = sum(delta,2)/N;
% covariance
Q = zeros(3);
for i=1:N
    Q = Q + (delta(:,i)-m)*(delta(:,i)-m)';
end
Q = Q/N-1;

figure;
plot3(delta(1,:),delta(2,:),delta(3,:),'.b');
axis equal;
end
