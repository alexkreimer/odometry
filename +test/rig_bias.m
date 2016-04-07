function rig_bias(num_pts, z_max, z_min)

% this test reproduces bias in rotation estimation for a stereo rig if
% there is a calibration error

close all;
dbstop if error;

dist_thr = 200;
num_frames = 2;

% KITTI intrinsics
K = [718.8560,  0,   607.1928;
    0,     718.8560, 185.2157;
    0,       0,      1.0000];

% distance between rig cameras
base = .5;

% right camera center as seen from the left camera
O2 = [base 0 0]';

T0 = [eye(3) O2; 0 0 0 1];

% left camera
P1 = K*[eye(3) zeros(3,1)];

% create our world
x_max = .02*z_max;
X(1,:) = -x_max + (x_max-(-x_max))*rand(1, num_pts); % 3d points, as seen in the world frame
X(2,:) = -x_max + (x_max-(-x_max))*rand(1, num_pts);
X(3,:) = z_min + (z_max-z_min)*rand(1, num_pts);

% allocate image projections
x = nan(2*num_frames, 2*num_pts);

% initial camera pose
T1 = [eye(3) zeros(3,1); 0 0 0 1];

% number of experiments
N = 10;

% preallocate epipoles
e      = nan(2,N,2);
e_gt = nan(2,num_frames-1);

% generate camera motion
R = rotx(.2*rand)*roty(.2*rand)*rotz(.2*rand);
t = [.5*rand .5*rand .5+rand]';

% orientation of the original frame in the new one
R_gt   = R';
t_gt   = -R'*t;
F = K'\util.skew(t_gt)*R_gt/K;

dt = [R t; 0 0 0 1];

for frame=1:num_frames
    % observations in the new image pair
    [x((2*frame-1):2*frame, 1:num_pts), visible1]              = util.project(P1, X, T1);
    [x((2*frame-1):2*frame,(num_pts+1):(2*num_pts)), visible2] = util.project(P1, X, T1*T0);
    
    if frame>1
        e_gt(:,frame-1) = util.h2e(K*(-R'*t));
    end
    
    T1 = dt*T1;
    assert(all([visible1,visible2]),'some points are invisible in the new pose of a rig');
end

% number of different measurements
n = 3;
names = {'left','right','stereo'};
r = nan(9,N,n);

% points
x1l = x(1:2, 1:num_pts); x1r = x(1:2,(num_pts+1):(2*num_pts));
x2l = x(3:4, 1:num_pts); x2r = x(3:4,(num_pts+1):(2*num_pts));
for i=1:N
    fprintf('experiment %d\n', i);
    
    if i>1 && 0
        x1l = x1l + 2*randn(size(x1l));
        x2l = x2l + 2*randn(size(x2l));
    end
    
    depths = X(3,:);
    % left
    tic;[T1, ~] = estimation.rel_motion_H(K,x1l,x2l,depths,base,'absRotInit',true,'F',F,'depth_thr',200);toc;
    r(:,i,1) = reshape(T1(1:3,1:3),[9 1]);

    % right
    tic;[T2, ~] = estimation.rel_motion_H(K,x1r,x2r,depths,base,'absRotInit',true,'F',F,'depth_thr',200);toc;
    r(:,i,2) = reshape(T2(1:3,1:3),[9 1]);
    
    % stereo
    xsl = [x1l x1r];
    xsr = [x2l x2r];
    tic;[T3, ~] = estimation.rel_motion_H(K,xsl,xsr,[depths, depths],base,'absRotInit',true,'F',F,'depth_thr',300);toc;
    r(:,i,3) = reshape(T3(1:3,1:3),[9 1]);
end

% compute rotation errors
dr = nan(4,N,n);
for i=1:N
    for j=1:n
        R1 = reshape(r(:,i,j),[3,3]);
        R2 = dt(1:3,1:3);
        dR = R1*R2';
        dr(:,i,j) = vrrotmat2vec(dR);
    end
end

figure; 
subplot(121); hold on;
title('rotation errors: [x y]');
for j=1:n
    scatter(dr(1,:,j),dr(2,:,j),'DisplayName',names{j});
end
legend show;

subplot(122); hold on;
title('rotation errors: \theta');
for j=1:n
    plot(dr(4,:,j),'o','DisplayName',names{j});
end
legend show;
