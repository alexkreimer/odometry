function chierality()
dbstop if error;
close all;

K = [718.8560,  0,      607.1928;
       0,     718.8560, 185.2157;
       0,       0,        1.0000];       % camera intrinsics

X = util.e2h(nan(3,100));

X(1,:) = 0;
X(2,:) = 0;
X(3,:) = 100*randn(1,100);

I = eye(3);
R = eye(3);
t = [1 0 0]';

P1 = K*[I zeros(3,1)];
P2 = K*[R' -R'*t];

x1 = P1*X;
x2 = P2*X;

detM1 = det(P1(1:3,1:3));
detM2 = det(P2(1:3,1:3));

visible = nan(1,length(X));
for i=1:length(X)
    visible(i) = X(4,i)*x1(3,i)*detM1 > 0 && X(4,i)*x2(3,i)*detM2>0;
end
figure; hold on;
plot(0,0,'or',t(1),0,'or');
plot(X(1,:),X(3,:));
figure; hold on;
plot(visible);
plot(X(3,:)>0,'or');
end
