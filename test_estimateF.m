function test_estimateF(num_pts, z_max, z_min)

OUT_DIR = '/home/kreimer/prj/odometry/debug';

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

X(1,:) = -3 + (3-(-3))*rand(1, num_pts); % 3d points, as seen in the world frame
X(2,:) = -3 + (3-(-3))*rand(1, num_pts);
X(3,:) = z_min + (z_max-z_min)*rand(1, num_pts);

x = nan(2*num_frames, 2*num_pts);     % image projections

T = [Id zeros(3,1); 0 0 0 1];
dt = [Id [0;0;1]; 0 0 0 1];
for frame=1:num_frames
    % observations in the new image pair
    
    [x((2*frame-1):2*frame, 1:num_pts), visible1]              = project(P1, X, T);
    [x((2*frame-1):2*frame,(num_pts+1):(2*num_pts)), visible2] = project(P1, X, T*T0);
    
    T = dt*T;
    assert(all([visible1,visible2]),'some points are invisible in the new pose of a rig');
end

% number of experiments
N = 50;

% preallocate
e = nan(2,N,2);


figure;
hold on;
title('feature tracks (both original and noisy)');

err   = nan(4,N);
r_err = nan(4,N);
t_err = nan(4,N);

for i=1:N
    fprintf('experiment %d\n', i);
    x1 = x(1:2, 1:num_pts);
    x2 = x(3:4, 1:num_pts);
    
    plot([x1(1,:); x2(1,:)], [x1(2,:); x2(2,:)]);
    
    if i>1
        x1 = x1 + randn(size(x1));
        x2 = x2 + randn(size(x2));
    end
    
    plot([x1(1,:); x2(1,:)], [x1(2,:); x2(2,:)]);
    
    [F, T, ~, ~,F_lin] = estimateF(x2, x1, K, 2, 'all');
    
    F = cat(3,F_lin, F);
    T = cat(3,[decompose_essential(K'*F_lin*K, K, [x1;x2]); 0 0 0 1], T);
    for j=1:4
        val = sampsonF(F(:,:,j), e2h(x1), e2h(x2));
        err(j,i) = sqrt(val'*val/N);
        T(1:3,4,j) = T(1:3,4,j)/norm(T(1:3,4,j));
        if T(3,4,j)>0
            T(:,:,j) = inv(T(:,:,j));
        end
        pose_error = T(:,:,j)*dt;
        r_err(j,i) = rot_error(pose_error);
        t_err(j,i) = trans_error(pose_error);
    end
end

figure; hold on; title('sampson error');
plot(err(1,:));
plot(err(2,:));
plot(err(3,:));
plot(err(4,:));
legend('8-pt', 'direct', 'vtheta', 'tangent');

figure; hold on;
plot(t_err(1,:));
plot(t_err(2,:));
plot(t_err(3,:));
plot(t_err(4,:));
title('translation error');
legend('8-pt', 'direct', 'vtheta', 'tangent');

figure; hold on;
plot(r_err(1,:));
plot(r_err(2,:));
plot(r_err(3,:));
plot(r_err(4,:));
title('rotation error');
legend('8-pt', 'direct', 'vtheta', 'tangent');