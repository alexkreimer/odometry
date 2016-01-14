function v = simple_H_inf()

K = [1 0 100; 0 1 100; 0 0 1];
K = [718.8560,  0,      607.1928;
       0,     718.8560, 185.2157;
       0,       0,        1.0000];       % camera intrinsics
X = randn(3,100);

I = eye(3);
R = rotx(rand*.2)*roty(rand*.2);
t = [.5*rand .5*rand .5+rand]';

x1 = util.h2e(K*[I zeros(3,1)]*util.e2h(X));
x2 = util.h2e(K*[R' -R'*t]*util.e2h(X));

H_inf = K*R'/K;

x1t_H = util.h2e(H_inf*util.e2h(x1));

delta = x2-x1t_H;

pd_x = fitdist(delta(1,:)','Normal');
pd_y = fitdist(delta(2,:)','Normal');

x_values = min(delta(1,:)'):max(delta(1,:)');
y_values = pdf(pd_x,x_values);
plot(x_values,y_values);

hold on;

x_values = min(delta(2,:)'):max(delta(2,:)');
y_values = pdf(pd_y,x_values);
plot(x_values,y_values);

end
