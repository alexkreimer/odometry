function estimate_motion1()

close all
dbstop if error;

KITTI_HOME = '/home/kreimer/KITTI/dataset';
KITTI_HOME = fullfile('F:', 'KITTI' , 'dataset');
DBG_DIR = fullfile('F:', 'debug');

image_dir  = fullfile(KITTI_HOME, 'sequences', '00');
poses_file = fullfile(KITTI_HOME, 'poses','00.txt');

% setup camera parameters (KITTI)
[P0, P1] = kitti_read_calib(image_dir);
poses_gt = kitti_read_poses(poses_file);
param = kitti_params(P0, P1);

gt = process_gt(poses_gt);
info = [];
% load features
load tracks;
num_frames = length(info);
poses1 = nan(4, 4, num_frames);
poses1(:,:,1) = inv([eye(3) zeros(3,1); 0 0 0 1]);
poses2 = poses1;
for i = 2:num_frames
    fprintf('processing frame %d\n', i);
    %[i1, i2] = read_kitti_images(image_dir, i);
    % collect features for the estimation
    cur = i; prv = i-1;
    c1 = info(cur).c1(:, info(cur).mt(1, :)); % current left
    c1p= info(prv).c1(:, info(cur).mt(2, :)); % previous left
    c2p= info(prv).c2(:, info(cur).mt(3, :)); % previous right
    
    mt = info(cur).mt;
    m12 = info(cur).m12;
    
    mt(4, :) = nan;
    for j = 1:length(mt)
        ind = find(m12(1, :) == mt(1, j));
        if ~isempty(ind)
            mt(4, j) = m12(2, ind);
        end
    end
    
    valid = ~isnan(mt(4, :));
    c1_ss = info(cur).c1(:, mt(1, valid));
    c2_ss = info(cur).c2(:, mt(4, valid));
    c1p_ss = info(prv).c1(:, mt(2, valid)); % previous left
    c2p_ss = info(prv).c2(:, mt(3, valid)); % previous right
    
    num_pts = length(c1_ss);
    X = nan(4, num_pts);
    for j = 1:num_pts
        X(:, j) = vgg_X_from_xP_nonlin([c1p_ss(:, j) c2p_ss(:,j)], {param.P1, param.P2});
    end
    X = h2e(X);
    
    [a_ss, ~, ~, ~, ~] = ransac_minimize_reproj(X, [c1_ss; c2_ss], param);
    T_ss = tr2mat(a_ss);
    E_ss = skew(T_ss(1:3,4))*T_ss(1:3,1:3);
    F_ss = inv(param.K')*E_ss*inv(param.K);
    
    poses2(:,:,i) = poses2(:,:,i-1)*inv(tr2mat(a_ss));
    
    % collect params
    params = struct('c1', c1,...
        'c1p', c1p,...
        'c2p', c2p,...
        't0', [param.base,0,0]',...   % stereo baseline
        'K', param.K);
    
    [a_est1, pout] = estimate_stereo_motion1(params);
    
    poses1(:,:,i) = poses1(:,:,i-1)*tr2mat(a_est1);
    
    % reconstruct F back from R,t
    T_rec = pout.est1.T;
    E_rec = skew(T_rec(1:3,4))*T_rec(1:3,1:3);
    F_rec = inv(params.K')*E_rec*inv(params.K);
    n = 2;
    
    pin = struct(...
        'x1', pout.est1.x1,...
        'x2', pout.est1.x2,...
        'inliers', pout.est1.inliers,...
        'i1', info(prv).i1,...
        'i2', info(cur).i1,...
        'i1_name', 'prev left',...
        'i2_name', 'cur left',...
        'ind1', prv,...
        'ind2', cur,...
        'DBG_DIR', DBG_DIR,...
        'K', params.K,...
        'dbg_save', 0);
    
    T_gt = gt(:, :, i);
    t_gt = T_gt(1:3, 4);
    R_gt = T_gt(1:3, 1:3);
    E_gt = skew(t_gt)*R_gt;
    F_gt = inv(params.K')*E_gt*inv(params.K);
    
    pin.F = {pout.est1.F, F_rec, F_gt, F_ss};
    pin.E = {pout.est1.E, E_rec, E_gt, E_ss};
    pin.T = {pout.est1.T, T_rec, T_gt, T_ss};
    
    pin.name = {'epipolar estimation', 'after reconstruction', 'gt', 'ss'};
    
    test_est_F(pin);
    
    pin = struct(...
        'x1', pout.est2.x1,...
        'x2', pout.est2.x2,...
        'inliers', pout.est2.inliers,...
        'i1', info(prv).i2,...
        'i2', info(cur).i1,...
        'i1_name', 'prev right',...
        'i2_name', 'cur left',...
        'ind1', prv,...
        'ind2', cur,...
        'DBG_DIR', DBG_DIR,...
        'K', params.K,...
        'dbg_save', 0);
    
    % reconstruct F back from R,t
    T_rec = [pout.est2.T; 0 0 0 1];
    E_rec = skew(T_rec(1:3,4))*T_rec(1:3,1:3);
    F_rec = inv(params.K')*E_rec*inv(params.K);
   
    T_gt(1,4) = T_gt(1,4) + param.base;
    t_gt = T_gt(1:3, 4);
    R_gt = T_gt(1:3, 1:3);
    E_gt = skew(t_gt)*R_gt;
    F_gt = inv(params.K')*E_gt*inv(params.K);
    
    T_ss(1,4) = T_ss(1,4) + param.base;
    E_ss = skew(T_ss(1:3,4))*T_ss(1:3,1:3);
    F_ss = inv(param.K')*E_ss*inv(param.K);
    
    pin.F = {pout.est2.F, F_rec, F_gt, F_ss};
    pin.E = {pout.est2.E, E_rec, E_gt, E_ss};
    pin.T = {pout.est2.T, T_rec, T_gt, T_ss};
    pin.name = {'epipolar estimation 2', 'after reconstruction', 'gt2', 'ss2'};
    
    test_est_F(pin);
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

