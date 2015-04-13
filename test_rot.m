close all;
t = pi/2;

ct = cos(t); st = sin(t);
R = [ct -st; st ct];
v = [0 1]';
figure;
plot(v(1),v(2),'xr');
rv = R*v;
hold on;
plot(rv(1),rv(2),'xb');
