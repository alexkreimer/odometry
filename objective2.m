function val = objective2(w, t0, x, ratio, sigma, h0, c)

% c is a 12 parameter vector that that parametrizes the motion of the
% camera in current frame (first 6 params) and previous frame (last 6
% params)

% t0 - baseline vector
% 
% x1 left/prev right feature locations for current frame
% x2 left/prev right feature locations for previous frame
% x3 left/left for current frame
% x4 left/left for previous frame

[q1, t1] = deal(c(1:3)', c(4:6)');
[q2, t2] = deal(c(7:9)', c(10:12)');

R1 = quaternion(param2quaternion(q1, h0(1,:)')).RotationMatrix;
R2 = quaternion(param2quaternion(q2, h0(2,:)')).RotationMatrix;

[x1, x2, x3, x4] = deal(x{1}, x{2}, x{3}, x{4});
v1 = sampson_err(R1, -t0 + t1, x1(4:6,:), x1(1:3,:));
v2 = sampson_err(R2, -t0 + t2, x2(1:3,:), x2(4:6,:));

v3 = sampson_err(R1, t1, x3(1:3,:), x3(4:6,:));
v4 = sampson_err(R2, t2, x4(1:3,:), x4(4:6,:));

delta = (abs(ratio)-(1+abs(norm(t1)/norm(t2))));

val = v1'*v1/length(v1) + ...
      v2'*v2/length(v2) + ...
      v3'*v3/length(v3) + ...
      v4'*v4/length(v4) + ...
      w*delta*delta/(sigma*sigma);

end