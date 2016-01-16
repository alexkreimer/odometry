function rel_motion(num_pts, z_max, z_min)

close all;
dbstop if error;

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
brk = 20;
num_pts = num_pts - brk;
X(1,:) = -x_max + (x_max-(-x_max))*rand(1, num_pts); % 3d points, as seen in the world frame
X(2,:) = -x_max + (x_max-(-x_max))*rand(1, num_pts);
X(3,:) = z_min + (z_max-z_min)*rand(1, num_pts);

X_close(1,:) = -x_max + (x_max-(-x_max))*rand(1, brk); % 3d points, as seen in the world frame
X_close(2,:) = -x_max + (x_max-(-x_max))*rand(1, brk);
X_close(3,:) = 2 + (5-2)*rand(1, brk);

X = [X X_close];

num_pts = num_pts + brk;

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
R_gt = R;
t_gt = t;
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

h = figure; subplot(231); hold on;
title('feature tracks (both original and noisy)');

names = {'F', 'FX', 'Fg', 'Hg', 'HX'};

r = nan(9,N,length(names));
e = nan(2,N,length(names));

for i=1:N
    fprintf('experiment %d\n', i);
    x1 = x(1:2, 1:num_pts); x1r = x(1:2,(num_pts+1):(2*num_pts));
    x2 = x(3:4, 1:num_pts); x2r = x(3:4,(num_pts+1):(2*num_pts));

    figure(h); subplot(231);
    for j=1:num_pts
        plot(x1(1,j), x1(2,j),'og');
        plot([x1(1,j); x2(1,j)], [x1(2,j); x2(2,j)],'r');
    end
    
    if i>1
        x1 = x1 + 1*randn(size(x1));
        x2 = x2 + 1*randn(size(x2));
    end
    
    plot([x1(1,:); x2(1,:)], [x1(2,:); x2(2,:)]);
    
    n = 1;
    % R from F decomposition, t null vector of F
    [F, T, ~, ~]  = estimation.estimateF(x1, x2, K, 2);
    R = T(1:3,1:3);
    e(:,i,n) = util.h2e(null(F));
    r(:,i,n) = reshape(R,[9 1]);

    n = n+1;
    % R as above, t by minimizing reprojection errors of 3d points
    t = estimation.trans_X(K,R,param.base,X,x2,x2r,-R_gt'*t_gt);
    e(:,i,n) = util.h2e(K*t);
    r(:,i,n) = reshape(R,[9 1]);
    
    n = n+1;
    % R as above, t by fitting epipolar lines
    H = K*R/K;
    t = estimation.trans_geom(K,H,x1,x2);
    e(:,i,n) = util.h2e(K*t);
    r(:,i,n) = reshape(R,[9 1]);    
    
    n = n+1;
    % R from H, t by fitting epipolar lines
    [T, inliers] = estimation.rel_motion_H(K, x1, x2, X(3,:), 1, 250, .01, e_gt);
    R = T(1:3,1:3);
    e(:,i,n) = util.h2e(K*T(1:3,4));
    r(:,i,n) = reshape(R,[9 1]);

    n = n+1;
    % R as above, t by minimizing reprojection errors of 3d points
    t = estimation.trans_X(K,R,param.base,X,x2,x2r,-R_gt'*t_gt);
    e(:,i,n) = util.h2e(K*t);
    r(:,i,n) = reshape(R,[9 1]);
end

figure(h);
% relative error in x and y
dist_fn = @(x1,x2) min([1, abs(x1-x2)/min(abs([x1,x2]))]);
rel_error = nan(2,N,2);
for i=1:N
    for j = 1:size(e,3)
        x1 = e_gt(1,1); x2 = e(1,i,j);
        y1 = e_gt(2,1); y2 = e(2,i,j);
        
        rel_error(1,i,j) = dist_fn(x1, x2);
        rel_error(2,i,j) = dist_fn(y1, y2);
    end
end

% relative error plot
subplot(232); hold on;
title('epipole worst relative coordinate error');
for j=1:size(e,3)
%    cur_rel_error = max(rel_error(:,:,j),[],1);
%    plot(cur_rel_error, 'DisplayName', sprintf('%s', names{j}));
    plot(rel_error(1,:,j)','DisplayName', sprintf('x-%s', names{j}));
%    plot(rel_error(2,:,j)','DisplayName', sprintf('y-%s', names{j}));
end
legend show;

% scatter the epipoles
subplot(233); hold on; title('epipoles');
for j = 1:size(e,3)
    scatter(e(1,:,j), e(2,:,j), 'DisplayName', names{j});
end
plot(e_gt(1),e_gt(2),'g*','DisplayName','gt');

legend show;

% fit normal distribution to the epipole estimations
subplot(234); hold on; title('epipole per-coordinate distributions');
for j= 1:size(e,3)
    data = e(:,:,j)';
    
    pd_x = fitdist(data(:,1), 'Normal');
%    pd_y = fitdist(data(:,2), 'Normal');
    
    x_values1 = min(data(:,1)):.1:max(data(:,1));
%    x_values2 = min(data(:,2)):.1:max(data(:,2));
    
    y1 = pdf(pd_x, x_values1);
%    y2 = pdf(pd_y, x_values2);
    plot(x_values1, y1,'LineWidth',2, 'DisplayName', sprintf('x-%s',names{j}));
%    plot(x_values2, y2,'LineWidth',2, 'DisplayName', sprintf('y-%s',names{j}));
    fprintf('x-distribution %s mean: %g, sigma: %g\n', names{j}, pd_x.mu, pd_x.sigma);
%    fprintf('y-distribution %s mean: %g, sigma: %g\n', names{j}, pd_y.mu, pd_y.sigma);
end
plot([e_gt(1,1);e_gt(1,1)],[0.002;0]);
%plot([e_gt(2,1);e_gt(2,1)],[0.002;0]);
legend show;

r_err = nan(2,N);
t_err = nan(2,N);

for i=1:N
    for j=1:size(e,3)
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
for j=1:size(e,3)
    plot(r_err(j,:),'DisplayName',names{j});
end
legend show;

subplot(236); hold on;
title({'average translation error', mean(t_err,2)});
for j=1:size(e,3)
    plot(t_err(j,:),'DisplayName',names{j});
end
legend show;
fprintf('real epipole %g %g\n', e_gt);

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