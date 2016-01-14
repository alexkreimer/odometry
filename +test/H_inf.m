function r_err_avg = H_inf(num_pts, z_max, z_min)

num_frames = 2;
K = [718.8560,  0,      607.1928;
    0,     718.8560, 185.2157;
    0,       0,        1.0000];       % camera intrinsics

% R0/t0 is the orientation/origin of the right camera as seen from the left camera
baseline            = [.5;0;0];
Id                  = eye(3);
T0                  = [Id, baseline; 0 0 0 1];
param.base          = norm(baseline);
param.calib.f       = K(1);
param.calib.cu      = K(1,3);
param.calib.cv      = K(2,3);
param.ransac_iter   = 1000;
param.model_size    = 3;
param.init          = true;
param.lm_max_iter   = 1000;
param.inlier_thresh = 3;
P1                  = K*[Id zeros(3,1)];     % left camera

x_max = .2*z_max;
X(1,:) = -x_max + (x_max-(-x_max))*rand(1, num_pts); % 3d points, as seen in the world frame
X(2,:) = -x_max + (x_max-(-x_max))*rand(1, num_pts);
X(3,:) = z_min + (z_max-z_min)*rand(1, num_pts);

% allocate image projections
x      = nan(2*num_frames, 2*num_pts);

% initial camera pose
T      = [Id zeros(3,1); 0 0 0 1];

% number of experiments
N      = 10;

% preallocate epipoles
e      = nan(2,N,2);
e_gt = nan(2,num_frames-1);

% generate camera motion
R      = rotx(.2*rand)*roty(.2*rand);
t      = [.5*rand .5*rand .5+rand]';
t      = t/norm(t);
%R = Id; t = [0 0 1]';
dt     = [R t; 0 0 0 1];

for frame=1:num_frames
    % observations in the new image pair
    [x((2*frame-1):2*frame, 1:num_pts), visible1]              = util.project(P1, X, T);
    [x((2*frame-1):2*frame,(num_pts+1):(2*num_pts)), visible2] = util.project(P1, X, T*T0);
    
    if frame>1
        e_gt(:,frame-1) = util.h2e(K*(-R'*t));
    end
    
    T = dt*T;
    assert(all([visible1,visible2]),'some points are invisible in the new pose of a rig');
end

h = figure; subplot(211); hold on;
title('feature tracks (both original and noisy)');

r_err = nan(1,N);

for i=1:N
    fprintf('experiment %d\n', i);
    x1 = x(1:2, 1:num_pts);
    x2 = x(3:4, 1:num_pts);

    figure(h); subplot(211);
    for j=1:num_pts
        plot(x1(1,j), x1(2,j),'og');
        plot([x1(1,j); x2(1,j)], [x1(2,j); x2(2,j)],'r');
    end
    
    if i>1
        x1 = x1 + 1*randn(size(x1));
        x2 = x2 + 1*randn(size(x2));
    end
    
    plot([x1(1,:); x2(1,:)], [x1(2,:); x2(2,:)]);
    
    R_est = estimation.H_inf_nonlin(K, x1, x2, X(3,:), 1, 250, .1);
    
    delta = R_est\R;
    pose_error = [delta zeros(3,1); 0 0 0 1];
    r_err(i) = util.rot_error(pose_error);
end

subplot(212); hold on; title('rotation error');
plot(r_err);
r_err_avg = mean(r_err);
