function teest_clsq()

N = 10;

x = rand(N,1);
w = rand(3,1);

x1 = x+.01*randn(N,1);
y1 = nan(N,1);
for i=1:N
    y1(i) = (-w(1)*x(i)-w(3))/w(2);
end

y1 = y1+.01*randn(N,1);

[c,n] = clsq([ones(N,1),x1,y1],2);

plotline(x1,y1,'o',c,n,'--');
