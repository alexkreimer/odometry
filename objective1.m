function val = objective1(t0, T1, T2, x1, x2, ratio, sigma, c1, c2)
% c1, c2 - initial values of the scale parameters (default =1)

% t0 - baseline vector
% T1, T2 sucessive left camera motions
% x1 left/prev right feature locations for T1
% x2 left/prev right feature locations for T2

% left camera translation vectors
q1 = T1(1:3,4);
q2 = T2(1:3,4);

% left camera rotations
R1 = T1(1:3,1:3);
R2 = T2(1:3,1:3);

% translations right camera vs left as seen from the right camera
t1 = -t0+c1*q1;
t2 = -t0+c2*q2;

v1 = sampson_err(R1,t1, x1(4:6,:), x1(1:3,:));
v2 = sampson_err(R2,t2, x2(1:3,:), x2(4:6,:));

delta = (abs(ratio)-(1+abs(c1/c2)));

val = v1'*v1/length(v1) + v2'*v2/length(v2) + delta*delta/(sigma*sigma);
end