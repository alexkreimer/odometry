function estimate_motion()

close all
load tracklets

KITTI_HOME = '/home/kreimer/KITTI/dataset';
KITTI_HOME = fullfile('F:', 'KITTI' , 'dataset');
DBG_DIR = fullfile('F:', 'debug');

image_dir  = fullfile(KITTI_HOME, 'sequences', '00');
poses_file = fullfile(KITTI_HOME, 'poses','00.txt');

% setup camera parameters (KITTI)
[P0, P1] = kitti_read_calib(image_dir);
poses_gt = kitti_read_poses(poses_file);


% baseline, focal, principal point
param.base = -P1(1,4)/P1(1,1);
param.K = P0(1:3,1:3);
param.calib.f = P0(1,1);
param.calib.cu = P0(1,3);
param.calib.cv = P0(2,3);
% this commands minimization procedure to first solve 3d-3d rigid motion
% and use it as an initial point for Gauss-Newton reprojection minimization
% procedure
param.init = true;
% fundamental
param.F = vgg_F_from_P(P0,P1);

param.feat_num = 4000;                       % number of corners to extract per image
param.feat_method = 'MinimumEigenvalue';     % corner extraction method (not used)
param.patchr = 3;                            % descriptors are of size (2*patchr+1)**2
param.searchr = 100;                         % search window radius for feature tracking

param.min_d = 4;                             % minimal acceptable disparity
param.min_disp = -inf;
param.ba_w = 3;                              % BA window

% reprojection error minimization error threshold
param.inlier_thresh = 1;
param.model_size = 3;
param.ransac_iter = 3000;
% observation dimension
param.obd = 4;
% camera parameter vector dimension
param.ad = 6;
% structure parameter vector dimension
param.bd = 3;
param.threshx = 100;
param.threshy = 2;
param.lm_max_iter = 100;

% process ground truth data, we want local transformations
step_num = size(poses_gt,3);
pgt = nan(3,4,step_num);
pgt(:,:,1) = poses_gt(:,:,1);
for i=2:step_num
    T1 = [poses_gt(:,:,i); 0 0 0 1];      % {i} -> {0}
    T2 = [poses_gt(:,:,i-1); 0 0 0 1];    % {i-1} -> {0}
    T1 = inv(T1);                         % {0} -> {i}
    T = T1*T2;                            % {i-1} -> {i}
    pgt(1:3,:,i) = T(1:3,:);
end

% baseline, focal, principal point
param.base = -P1(1,4)/P1(1,1);
param.K = P0(1:3,1:3);

max_step = max([tracklets1.step]);

poses1(:,:,1) = tr2mat(zeros(6,1));
poses2(:,:,1) = tr2mat(zeros(6,1));

[pi1,pi2] = read_kitti_images(image_dir, 1);

for i = 2:10
    fprintf('step %d\n', i);
    [i1,i2] = read_kitti_images(image_dir, i);
    [x, px, Xp] = tracklets1.get_matched(tracklets2, i);
    % plot_circles(px, x, i1, pi1, i2, pi2);
    num_pts = size(x,2);
    x = [x(1:2,:) x(3:4,:); px(1:2,:) px(3:4,:)];
    a_est = estimate_stereo_motion(x, param.K, num_pts, eye(3), [param.base,0,0]',...
                                   'DBG_DIR',DBG_DIR, 'gt', pgt(:,:,i), 'ind', i, 'i1', i1, 'pi1', pi1, 'pi2', pi2);
    T = tr2mat(a_est);
    poses1(:,:,i) = T*poses1(:,:,i-1);
    figure; scatter3(Xp(1,:), Xp(2,:), Xp(3,:));
    % estimate params
    [a_est, inliers, tr0, predict, rms] = ransac_minimize_reproj(Xp, [x(1:2,1:num_pts); x(3:4, (num_pts+1):end)], param);
    a_est = tinv(a_est);
    T = tr2mat(a_est);
    poses2(:,:,i) = T*poses2(:,:,i-1);
    
    pi1 = i1; pi2 = i2;
end

savePoses('00_f.txt', poses1);
savePoses('00_3d.txt', poses2);
end

function plot_circles(px, x, i1, pi1, i2, pi2)
% px = [px1;px2]
%  x = [ x1; x2]

px(2,:) = px(2,:) + size(i1,1);
 x(3,:) =  x(3,:) + size(i1,2);
px(3,:) = px(3,:) + size(i1,2);
px(4,:) = px(4,:) + size(i1,1);

im = [i1, i2; pi1, pi2];
figure; imshow(im); hold on; title('circular matches');
for i=1:size(x,2)
    plot([x(1,i) x(3,i) px(3,i) px(1,i) x(1,i)],...
         [x(2,i) x(4,i) px(4,i) px(2,i) x(2,i)]);
end
hold off;
end

