function estimate_motion1()

close all
dbstop if error;

KITTI_HOME = '/home/kreimer/KITTI/dataset';
KITTI_HOME = '/media/kreimer/my_drive/record_20150720/dataset/';
DBG_DIR = fullfile('/home/kreimer/tmp', 'debug');

sequence = 's00';
image_dir  = fullfile(KITTI_HOME, 'sequences', sequence);

% setup camera parameters (KITTI)
[P0, P1] = kitti_read_calib(image_dir);
param = kitti_params(P0, P1);

info = struct('i1',[],'i2',[],'c1',[],'c2',[],'f1',[],'f2',[],'m12',[],'m11p',[],'m22p',[],'m12p',[],'m21p',[],'mc',[],'mt',[]);
% load features
num_frames = 100;
for i=1:num_frames
    load(['tracks/tracks_', sequence, '_', int2str(i), '.mat']);
    info(i) = val;
end

poses1 = nan(4, 4, num_frames);
poses1(:,:,1) = inv([eye(3) zeros(3,1); 0 0 0 1]);
as = nan(6,num_frames);
for i = 2:num_frames
    fprintf('processing frame %d\n', i);
    [i1, i2] = read_kitti_images(image_dir, i);
    % collect features for the estimationb
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
    
    figure; hist(X(3,:));
    
    param.init = 0;
    [a_ss, inliers, ~, ~, ~] = ransac_minimize_reproj(X, [c1_ss; c2_ss], param);
    as(:,i) = a_ss;
    figure;
    imshow(i1);
    hold on;
    plot(c1_ss(1,inliers),c1_ss(2,inliers),'.r');
    T_ss = tr2mat(a_ss);
    E_ss = skew(T_ss(1:3,4))*T_ss(1:3,1:3);
    F_ss = inv(param.K')*E_ss*inv(param.K);
    poses1(:,:,i) = poses1(:,:,i-1)*tr2mat(a_ss);
end
savePoses([sequence, '_3d.txt'], poses1);
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