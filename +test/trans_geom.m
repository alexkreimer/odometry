function trans_geom()
close all;
dbstop if error;

K = [718.8560,  0,      607.1928;
       0,     718.8560, 185.2157;
       0,       0,        1.0000];       % camera intrinsics
   
N = 100;
X = util.e2h(10+10*randn(3,N));

I = eye(3);
R = rotx(randn*.2)*roty(randn*.2);
t = [.5*randn .5*randn 2+.5*randn]';

x1 = K*[I zeros(3,1)]*X;
x2 = K*[R' -R'*t]*X;

H = K*R'/K;

Hx1   = H*x1;

e = estimation.trans_geom(K,H,Hx1,x2);
end
