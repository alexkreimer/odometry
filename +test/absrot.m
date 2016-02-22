function absrot()

X = rand(3,10);

K = [718.8560,  0,      607.1928;
    0,     718.8560, 185.2157;
    0,       0,        1.0000];

R = rotx(.2*rand)*roty(.2*rand);

x1 = util.h2e(K*X);
x2 = util.h2e(K*R*X);

p1 = K\util.e2h(x1);
p2 = K\util.e2h(x2);


for i=1:length(p1)
    p1(:,i) = p1(:,i)/norm(p1(:,i));
    p2(:,i) = p2(:,i)/norm(p2(:,i));
end

[regParams,~,~] = util.absor(p1,p2,'doTrans',false);

p2-regParams.R*p1