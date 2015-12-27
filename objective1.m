function val = objective1(w, K, t0, T, x, ratio, sigma, c)
% c1, c2 - initial values of the scale parameters (default =1)

% t0 - baseline vector
% T1, T2 sucessive left camera motions
% x1 left/prev right feature locations for T1
% x2 left/prev right feature locations for T2

[x1, x2] = deal(x{1}, x{2});
[T1, T2] = deal(T{1}, T{2});

% left camera translation vectors
q1 = T1(1:3,4);
q2 = T2(1:3,4);

% left camera rotations
R1 = T1(1:3,1:3);
R2 = T2(1:3,1:3);

% translations right camera vs left as seen from the right camera
t1 = -t0+c(1)*q1;
t2 = -t0+c(2)*q2;

t1 = q1;
t2 = q2;

F1 = K'\(skew(t1)*R1)/K;
F2 = K'\(skew(t2)*R2)/K;

v1 = sampsonF(F1, x1(1:3,:), x1(4:6,:));
v2 = sampsonF(F2, x2(1:3,:), x2(4:6,:));

delta = (abs(ratio)-(1+abs(c(1)/c(2))));

val = [v1/sqrt(length(v1)); v2/sqrt(length(v2)); w*delta/(sigma*sigma)];
end