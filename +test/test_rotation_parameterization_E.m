function test_rotation_parameterization_E(num_pts, z_max, z_min)

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
figure; 
T = [Id zeros(3,1); 0 0 0 1];
dt = [Id [0;0;1]; 0 0 0 1];
Q      = rotx(.2*rand)*roty(.2*rand);
t      = [.5*rand .5*rand .5+rand]';
dt     = [Q t; 0 0 0 1];
for frame=1:num_frames
    % observations in the new image pair
    
    [x((2*frame-1):2*frame, 1:num_pts), visible1]              = util.project(P1, X, T);
    [x((2*frame-1):2*frame,(num_pts+1):(2*num_pts)), visible2] = util.project(P1, X, T*T0);
    
    T = dt*T;
    assert(all([visible1,visible2]),'some points are invisible in the new pose of a rig');
end

% number of experiments
N = 10;

% preallocate
e = nan(2,N,2);
m=2; n=2;
subplot(m,n,4);
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
    
    [F, T, ~, inliers, residual, F_lin] = estimation.estimateF(x2, x1, K, 2, 'all');
    
    F = cat(3,F_lin,F);
    T = cat(3,[util.decompose_essential(K'*F_lin*K, K, [x1;x2]); 0 0 0 1], T);
    residual = [estimation.sampsonF(F_lin, util.e2h(x1(:,inliers')), util.e2h(x2(:,inliers'))) residual];
        
    for j=1:size(F,3)
        err(j,i) = sqrt(residual(:,j)'*residual(:,j)/size(residual,2));
        T(1:3,4,j) = T(1:3,4,j)/norm(T(1:3,4,j));
        if T(3,4,j)>0
            T(:,:,j) = inv(T(:,:,j));
        end
        pose_error = T(:,:,j)*dt;
        r_err(j,i) = mod(util.rot_error(pose_error),pi);
        r_err(j,i) = min(pi-r_err(j,i),r_err(j,i));
        t_err(j,i) = util.trans_error(pose_error);
    end
end

names = {'8-pt', 'direct', 'v-theta', 'tangent'};
subplot(m,n,1);
hold on;
avg_err = sum(err,2)/size(err,2);
title(sprintf('sampson error: %.2g %.2g %.2g %.2g',avg_err));
for i = 1:size(err,1)
    plot(err(i,:), 'DisplayName', names{i});
end
legend show;

subplot(m,n,2);
hold on;
for i = 1:size(t_err,1)
    plot(t_err(i,:),'DisplayName', names{i});
end
title(sprintf('t-err: %.2g %.2g %.2g %.2g',sum(t_err,2)/size(t_err,2)));
legend show;

subplot(m,n,3);
hold on;
for i = 1:size(r_err,1)
    plot(r_err(i,:),'DisplayName', names{i});
end
title(sprintf('r-err: %.2g %.2g %.2g %.2g',sum(r_err,2)/size(r_err,2)));

legend show;