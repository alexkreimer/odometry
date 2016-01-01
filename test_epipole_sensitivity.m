function test_epipole_sensitivity(num_pts, z_max, z_min)

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

x_max = .2*z_max;
X(1,:) = -x_max + (x_max-(-x_max))*rand(1, num_pts); % 3d points, as seen in the world frame
X(2,:) = -x_max + (x_max-(-x_max))*rand(1, num_pts);
X(3,:) = z_min + (z_max-z_min)*rand(1, num_pts);

% allocate image projections
x      = nan(2*num_frames, 2*num_pts);     

% initial camera pose
T      = [Id zeros(3,1); 0 0 0 1];

% number of experiments
N      = 50;

% preallocate epipoles
e      = nan(2,N,2);
e_true = nan(2,num_frames);

% generate camera motion
Q      = rotx(.2*rand)*roty(.2*rand);
t      = [.5*rand .5*rand .5+rand]';
%Q = Id; t = [-.1 0 1]';

dt     = [Q t; 0 0 0 1];

for frame=1:num_frames
    % observations in the new image pair
    [x((2*frame-1):2*frame, 1:num_pts), visible1]              = project(P1, X, T);
    [x((2*frame-1):2*frame,(num_pts+1):(2*num_pts)), visible2] = project(P1, X, T*T0);

    if frame>1
        e_true(:,frame-1) = project(P1, T(1:3,4));
    end

    T = dt*T;
    assert(all([visible1,visible2]),'some points are invisible in the new pose of a rig');
end

figure; subplot(231); hold on;
title('feature tracks (both original and noisy)');

r = nan(9,N,2);
for i=1:N
    fprintf('experiment %d\n', i);
    x1 = x(1:2, 1:num_pts);
    x2 = x(3:4, 1:num_pts);
    
    for j=1:num_pts
        plot(x1(1,j), x1(2,j),'og');
        plot([x1(1,j); x2(1,j)], [x1(2,j); x2(2,j)],'r');
    end
    
    if i>1
        x1 = x1 + 2*randn(size(x1));
        x2 = x2 + 2*randn(size(x2));
    end
    
    plot([x1(1,:); x2(1,:)], [x1(2,:); x2(2,:)]);

    [F, T, ~, ~] = estimateF(x1, x2, K, 2);
    [e(:,i,2), R] = epipole_from_H(K, x1, x2, X(3,:), 1);
    e(:, i, 1) = h2e(null(F));
    
    r(:,i,1) = reshape(T(1:3,1:3), [9 1]);
    r(:,i,2) = R(:);
end

% relative error in x and y
dist_fn = @(x1,x2) min([1, abs(x1-x2)/min(abs([x1,x2]))]);
rel_error = nan(2,N,2);
for i=1:N
    for j = 1:2
        x1 = e_true(1,1); x2 = e(1,i,j);
        y1 = e_true(2,1); y2 = e(2,i,j);
        
        rel_error(1,i,j) = dist_fn(x1, x2);
        rel_error(2,i,j) = dist_fn(y1, y2);
    end
end

% relative error plot
subplot(232); hold on;
title('epipole relative coordinate error');
plot(rel_error(:,:,1)');
plot(rel_error(:,:,2)');
legend('x-F','y-F','x-H','y-H');

% scatter the epipoles
subplot(233); hold on;
plot(e(1,:,1), e(2,:,1), '.r');
plot(e(1,:,2), e(2,:,2), 'ob');
plot(e_true(1,1), e_true(2,1), '*g');
legend('epipole(s) from F', 'epipole(s) from H', 'real epipole');

% fit normal distribution to the epipole estimations
data1 = e(:,:,1)';
data2 = e(:,:,2)';

pd1 = fitdist(data1(:,1), 'Normal');
pd2 = fitdist(data1(:,2), 'Normal');
pd3 = fitdist(data2(:,1), 'Normal');
pd4 = fitdist(data2(:,2), 'Normal');

x_values1 = min(data1(:,1)):.1:max(data1(:,1));
x_values2 = min(data1(:,2)):.1:max(data1(:,2));
x_values3 = min(data2(:,1)):.1:max(data2(:,1));
x_values4 = min(data2(:,2)):.1:max(data2(:,2));

y1 = pdf(pd1, x_values1);
y2 = pdf(pd2, x_values2);
y3 = pdf(pd3, x_values3);
y4 = pdf(pd4, x_values4);

subplot(234); hold on; title('epipole per-coordinate distributions');
plot(x_values1, y1,'LineWidth',2);
plot(x_values2, y2,'LineWidth',2);
plot(x_values3, y3,'LineWidth',2);
plot(x_values4, y4,'LineWidth',2);
legend('x of epiple from F', 'y of epipole from F', 'x of epipole from H', 'y of epipole from H');

r_err = nan(2,N);
t_err = nan(2,N);

for i=1:N
    r1 = reshape(r(:,i,1),[3,3]); t1 = K\e2h(e(:,i,1)); t1 = t1/norm(t1);
    T1 = [r1 t1; 0 0 0 1];
    r2 = reshape(r(:,i,2),[3,3]); t2 = K\e2h(e(:,i,2)); t2 = t2/norm(t2);
    T2 = [r2 t2; 0 0 0 1];
    
    pose_error = T1\dt;
    r_err(1,i) = mod(rot_error(pose_error), pi);
    r_err(1,i) = min(r_err(1,i),pi-r_err(1,i));
    t_err(1,i) = trans_error(pose_error);
    
    pose_error = T2\dt;
    r_err(2,i) = rot_error(pose_error);
    t_err(2,i) = trans_error(pose_error);
end

subplot(235); hold on; title('rotation error');
plot(r_err(1,:));
plot(r_err(2,:));
legend('F','H');

subplot(236); hold on; title('translation error');
plot(t_err(1,:));
plot(t_err(2,:));
legend('F','H');

fprintf('x-distribution-F mean: %g, sigma: %g\n', pd1.mu, pd1.sigma);
fprintf('y-distribution-F mean: %g, sigma: %g\n', pd2.mu, pd2.sigma);
fprintf('x-distribution-H mean: %g, sigma: %g\n', pd3.mu, pd3.sigma);
fprintf('y-distribution-H mean: %g, sigma: %g\n', pd4.mu, pd4.sigma);

fprintf('real epipole %g %g\n', e_true);

end

function val = sampson_error(E, x1, x2)
% x1, x2 inlier coordinates, s.t. x2'Ex1=0

N = size(x1,2);

val = nan(N,1);
for i=1:N
    % current point
    p1 = x1(:,i);
    p2 = x2(:,i);
    
    if length(p1) == 2
        p1 = [p1;1];
    end
    
    if length(p2) == 2
        p2 = [p2;1];
    end

    p2tE = p2'*E;
    Ep1  = E*p1;
    
    val(i) = p2'*E*p1/sqrt(p2tE(1)^2 + (p2tE(1))^2 + (Ep1(1))^2 + (Ep1(2))^2);
end
end