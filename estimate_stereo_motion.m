% This function estimates the motion of a stereo rig between 2 pairs of stereo frames

function a_est = estimate_stereo_motion(x,K,X,num_pts,R0,t0)
% x(1:2,1:num_pts) projections in the left camera, first pose
% x(1:2,num_pts+1:2*num_pts) projections in the right camera, first pose
% x(3:4,1:num_pts) projections if the left camera, 2nd pose
% x(3:4,num_pts+1:2*num_pts) projections in the right camera, 2nd pose
% K intrinsic calibration matrix
% 3d point(s) visible in the first camera, used to resolve the twisted pair
% ambiguity
% R0,t0 extrinsics of the stereo rig

% Estimate the motion of the left camera.  We determine all motion parametes
% except the singed magnitue of the translation
x1 = [x(1:2,1:num_pts); ones(1,num_pts)];
x2 = [x(3:4,1:num_pts); ones(1,num_pts)];
F1 = fundmatrix(x1,x2); % x2'*F*x1
E1 = K'*F1*K;
s1 = invE(E1,X(:,1),K);

t1_est = s1(1:3,4);
t1_est = t1_est/norm(t1_est);
if t1_est(3)<0
    % we assume that the camera is moving forward
    t1_est = -t1_est;
end
R1_est = s1(1:3,1:3);

% Estimate motion between initial position of the right camera and current
% position of the left camera
x1 = [x(1:2,(num_pts+1):(2*num_pts)); ones(1,num_pts)];
x2 = [x(3:4,1:num_pts); ones(1,num_pts)];
F2 = fundmatrix(x1,x2); % x2'*F*x1
E2 = K'*F2*K;
s2 = invE(E2,X(:,1),K);
t2_est = s2(1:3,4);
t2_est = t2_est/norm(t2_est);
if t2_est(3)<0
    % going forward
    t2_est = -t2_est;
end

t2_est = R0'*t2_est;     % t2 in the frame of the left camera

% now solve for scale of the motion
A = [t1_est, -t2_est];
b = t0;
[U,S,V] = svd(A);
b = U'*b;
d = diag(S);
y = nan(size(A,2),1);
for j=1:length(d)
    y(j) = b(j)/d(j);
end
c = V*y;

% form the output
t1_est = t1_est*c(1);
[rx,ry,rz] = decompose_rotation(R1_est);
a_est =  [rx,ry,rz,t1_est'];
end
