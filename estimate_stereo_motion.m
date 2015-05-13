% This function estimates the motion of a stereo rig between 2 pairs of stereo frames

function a_est = estimate_stereo_motion(x,K,num_pts,R0,t0,a0)
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

% Estimate 1st essential: x2'*F*x1
E1 = K'*fundmatrix(x1,x2)*K; 

% Decompose essential
s1 = decompose_essential(E1,K,x(1:4,1:num_pts));

% origin and orientation of the left_cam(1) as seen from left_cam(2)
t1 = s1(1:3,4);
R1 = s1(1:3,1:3);

% if we have ground truth, verify estimation
if nargin>5
    T1 = inv([s1; 0 0 0 1]);
    T  = tr2mat(a0);
    assertT(T1,T);
end

% Estimate motion between initial position of the right camera and current
% position of the left camera
x1 = x(1:2,(num_pts+1):(2*num_pts));
x2 = x(3:4,1:num_pts);

E2 = K'*fundmatrix(e2h(x1),e2h(x2))*K; % x2'*F*x1
s2 = decompose_essential(E2,K,[x1;x2]);

if nargin>5
    % pose of lcam2 in rcam1
    T2 = [s2; 0 0 0 1];
    
    % rcam1 described in lcam1
    T0 = [R0, t0;0 0 0 1];
    
    % lcam2 in rcam1
    T = T0\T;
    T = unit_t(T);
    
    assertT(T2,T,'inv');
end

t2 = s2(1:3,4);
R2 = s2(1:3,1:3);

t2 = -R2'*t2;     % t2 in the frame of the left camera
t1 = -R1'*t1;

% now solve for scale of the motion
A = [t1, -t2];
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
t1 = t1*c(1);
T1 = [R1' t1;0 0 0 1];

% if we have ground truth, verify estimation
if nargin>5
    % ground truth is a pose of cam2 in cam1
    assertT(T1,tr2mat(a0));
end

a_est = mat2tr(T1);

end

function T = unit_t(T)
T(1:3,4) = T(1:3,4)/norm(T(1:3,4));
end
