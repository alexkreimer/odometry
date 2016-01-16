function v = simple_trans_X()

K = [718.8560,  0,      607.1928;
       0,     718.8560, 185.2157;
       0,       0,        1.0000];       % camera intrinsics
X = randn(3,5);
base = [.5;0;0];

I = eye(3);
R = rotx(rand*.2)*roty(rand*.2);
t = [.5*rand .5*rand .5+rand]';

P1l = K*[I zeros(3,1)];
P1r = K*[I base];

x1l = util.h2e(P1l*util.e2h(X));
x1r = util.h2e(P1r*util.e2h(X));

X = util.triangulate_chieral(x1l,x1r,P1l,P1r);

x2l = util.h2e(K*[R' -R'*t]*util.e2h(X));
x2r = util.h2e(K*[R' -R'*(t+base)]*util.e2h(X));

x2l = x2l + randn(size(x2l));
x2r = x2r + randn(size(x2r));

t_est = estimation.trans_X(K,R',base(1),X,x2l,x2r,-R'*t);
norm(t_est-(-R'*t/norm(-R'*t)))

end
