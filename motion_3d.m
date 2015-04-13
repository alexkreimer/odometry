function motion_3d()
close all;

X = [0 0 0; 1 0 0]';

p = [.1,.1,.1,1,1,1];

T = tr2mat(p);
X1 = h2e(T*e2h(X));

A = [X,X1]';
rank(A)

figure;
scatter3(X(1,:),X(2,:),X(3,:));
hold on;
scatter3(X1(1,:),X1(2,:),X1(3,:));

