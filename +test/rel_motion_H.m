function rel_motion_H(num_pts, z_max, z_min, titl)

close all;
dbstop if error;

num_frames = 2;
K = [718.8560,  0,      607.1928;
    0,     718.8560, 185.2157;
    0,       0,        1.0000];       % camera intrinsics

K = [718.8560,  0,      0;
    0,     718.8560,    0;
    0,       0,        1.0000]; 
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

x_max = 10;
X(1,:) = x_max*rand(1, num_pts); % 3d points, as seen in the world frame
X(2,:) = x_max*rand(1, num_pts);
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
R      = rotx(.2*rand)*roty(.2*rand)*rotz(.2*rand);
t      = [.5*rand .5*rand .5+rand]';
R_gt   = R';
t_gt   = -R'*t;
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

h = figure; subplot(231); hold on;
title('feature tracks (both original and noisy)');

% number of different measurements
n = 2;
names = {'F','no-F'};

e = nan(2,N,n);
r = nan(9,N,n);

x1 = x(1:2, 1:num_pts); x1r = x(1:2,(num_pts+1):(2*num_pts));
x2 = x(3:4, 1:num_pts); x2r = x(3:4,(num_pts+1):(2*num_pts));

F = K'\util.skew(t_gt)*R_gt/K;
    
figure;
hold on;

e = util.h2e(null(F'));


H = K*R_gt/K;
Hx1 = util.h2e(H*util.e2h(x1));



epilines = epipolarLine(F,x1');
points = lineToBorderPoints(epilines, [2500 2500]);
line(points(:, [1,3])', points(:, [2,4])', 'LineWidth',2);
h1 = scatter(Hx1(1,:),Hx1(2,:), 40, 'filled');
h2 = scatter(x1(1,:),x1(2,:), 40, 'filled');
h3 = scatter(x2(1,:),x2(2,:), 40, 'filled');
h4 = plot(e(1),e(2),'r*');

legend([h2, h3, h4, h1], {'x','x''','Epipole','Hx'});

print('-depsc','-tiff','feature_motion_far')
return;

for i=1:N
    fprintf('experiment %d\n', i);
    
    figure(h); subplot(231);
    for j=1:num_pts
        plot(x1(1,j), x1(2,j),'og');
        plot([x1(1,j); x2(1,j)], [x1(2,j); x2(2,j)],'r');
    end
    
    if i>1
        x1 = x1 + 2*randn(size(x1));
        x2 = x2 + 2*randn(size(x2));
    end

    disp(names{1});
    save('R','R_gt');
    tic;[T, ~] = estimation.rel_motion_H(K,x1,x2,X(3,:),1,'depth_thr',100,'absRotInit',true);toc;
    R = T(1:3,1:3);
    e(:,i,1) = util.h2e(K*t_gt);
    r(:,i,1) = reshape(R,[9 1]);

    disp(names{2});
    tic;[T, ~] = estimation.rel_motion_H(K,x1,x2,X(3,:),1,'depth_thr',100,'absRotInit',false);toc;
    R = T(1:3,1:3);
    %e(:,i,2) = util.h2e(K*T(1:3,4));
    e(:,i,2) = util.h2e(K*t_gt);
    r(:,i,2) = reshape(R,[9 1]);
end

figure(h);
% relative error in x and y
dist_fn = @(x1,x2) min([1,abs(x1-x2)/min(abs([x1,x2]))]);
rel_error = nan(2,N,2);

x1 = e_gt(1,1);
y1 = e_gt(2,1);
for i=1:N
    for j = 1:n
        x2 = e(1,i,j);
        y2 = e(2,i,j);

        rel_error(1,i,j) = dist_fn(x1,x2);
        rel_error(2,i,j) = dist_fn(y1,y2);
    end
end

% relative error plot
subplot(232); hold on;
title('epipole worst relative error');
for j=1:n
    plot(rel_error(1,:,j)','DisplayName', sprintf('x-%s',names{j}), 'LineWidth', 5);
end
legend show;

% scatter the epipoles
subplot(233); hold on; title('epipoles');
for j=1:n
    scatter(e(1,:,j), e(2,:,j),'DisplayName',sprintf('estimates-%s',names{j}));
end
plot(e_gt(1),e_gt(2),'g*','DisplayName','gt');
legend show;

% fit normal distribution to the epipole estimations
subplot(234); hold on; title('epipole per-coordinate distributions');
for j=1:n
    data = e(:,:,j)';

    pd_x = fitdist(data(:,1), 'Normal');
    pd_y = fitdist(data(:,2), 'Normal');

    x_values1 = min(data(:,1)):.1:max(data(:,1));
    x_values2 = min(data(:,2)):.1:max(data(:,2));

    y1 = pdf(pd_x, x_values1);
    y2 = pdf(pd_y, x_values2);
    plot(x_values1, y1,'LineWidth',2, 'DisplayName', 'x');
    plot(x_values2, y2,'LineWidth',2, 'DisplayName', 'y');
    fprintf('x-distribution mean: %g, sigma: %g\n', pd_x.mu, pd_x.sigma);
    fprintf('y-distribution mean: %g, sigma: %g\n', pd_y.mu, pd_y.sigma);
end
plot([e_gt(1,1);e_gt(1,1)],[0.002;0]);
legend show;

r_err = nan(n,N);
t_err = nan(n,N);
for i=1:N
    for j=1:n
        R = reshape(r(:,i,j),[3,3]);
        t = K\util.e2h(e(:,i,j));
        t = t/norm(t);
        T = [R t; 0 0 0 1];

        pose_error = T\dt;
        r_err(j,i) = mod(util.rot_error(pose_error), pi);
        r_err(j,i) = min(r_err(j,i),pi-r_err(j,i));
        t_err(j,i) = util.trans_error(pose_error);
    end
end

subplot(235); hold on;
title({'average rotation errors:',mean(r_err,2)});
for j=1:n
    plot(r_err(j,:),'DisplayName',names{j});
end
ylim([0 1])
legend show;

subplot(236); hold on;
title({'average translation error', mean(t_err,2)});
for j=1:n
    plot(t_err(j,:), 'DisplayName',names{j});
end
ylim([0 1])
legend show;
fprintf('real epipole %g %g\n', e_gt);
end