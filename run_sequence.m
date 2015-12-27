function run_sequence(sequence, num_frames)

close all
dbstop if error;

DATA_ROOT  = '/media/kreimer/my_drive/KITTI/';
KITTI_HOME = fullfile(DATA_ROOT, 'dataset');
RESULT_DIR = fullfile(DATA_ROOT, 'results');
DBG_DIR    = fullfile(DATA_ROOT, 'debug');

image_dir  = fullfile(KITTI_HOME, 'sequences', sequence);
poses_file = fullfile(KITTI_HOME, 'poses', [sequence, '.txt']);

if nargin<2
    num_frames = 10;
end

% setup camera parameters (KITTI)
[P0, P1] = kitti_read_calib(image_dir);
poses_gt = kitti_read_poses(poses_file);
param = kitti_params(P0, P1);

gt = process_gt(poses_gt);
tracks = struct('i1',[],...
                'i2',[],...
                'c1',[],...
                'c2',[],...
                'f1',[],...
                'f2',[],...
                'm12',[],...
                'm11p',[],...
                'm22p',[],...
                'm12p',[],...
                'm21p',[],...
                'mc',[],...
                'mt',[]);

% load features
for i=1:num_frames
    load([DATA_ROOT, 'tracks/', sequence, '/frame_', num2str(i), '.mat']);
    tracks(i) = eval('info');
end

info = tracks;
poses1 = nan(4, 4, num_frames);
poses1(:,:,1) = inv([eye(3) zeros(3,1); 0 0 0 1]);
poses2 = poses1;

M = tracks_collect(info, 1);
stats = struct('support_size',[],'ratio',[],'sigma', []);

results = {'ss', 'no_opt'};

for i = 2:num_frames
    fprintf('processing frame %d\n', i);
    stats(i).support_size = size(info(i).mt,2);
    
    %[i1, i2] = read_kitti_images(image_dir, i);
    
    % collect features for the estimation
    cur = i;
    prv = i-1;
    c1  = info(cur).c1(:, info(cur).mt(1, :)); % current left
    c1p = info(prv).c1(:, info(cur).mt(2, :)); % previous left
    c2p = info(prv).c2(:, info(cur).mt(3, :)); % previous right

    % our algorithm only uses 1.5 stereo pair, while stereoscan needs all 4
    % images.  This functions completes the match into the current right
    % image
    mt  = info(cur).mt;
    m12 = info(cur).m12;    
    mt = complete_circle(mt, m12);

    % collect feature points for stereo-scan
    valid  = ~isnan(mt(4, :));
    c1_ss  = info(cur).c1(:, mt(1, valid));
    c2_ss  = info(cur).c2(:, mt(4, valid));
    c1p_ss = info(prv).c1(:, mt(2, valid)); % previous left
    c2p_ss = info(prv).c2(:, mt(3, valid)); % previous right
    % produce 3d
    X = triangulate_points(c1p_ss, c2p_ss, nnz(valid), param);
    % optimization
    [a_ss, ~, ~, ~, ~] = ransac_minimize_reproj(X, [c1_ss; c2_ss], param);
    % save the results
    T_ss = tr2mat(a_ss);
    E_ss = skew(T_ss(1:3,4))*T_ss(1:3,1:3);
    F_ss = inv(param.K')*E_ss*inv(param.K);
    stats(i).ss.T = tr2mat(a_ss);
    
    % collect params
    params = struct('c1',  c1, ...
                    'c1p', c1p,...
                    'c2p', c2p,...
                    'tracksx', [],...
                    'tracksy', [],...
                    't0', [param.base,0,0]',...   % stereo baseline
                    'K', param.K);

    % M is a cell-array; M{j} holds all tracklets of length j, found in
    % frame i. A = M{j} is an j\times N array (N is the number of tracklets)
    % Each column of A contains indices of features that comprise the
    % tracklet. A(k,l) is the index of the feature in frame i-j+k
    M = tracks_collect(info, i, M);
    
    pout = estimate_stereo_motion_new(params);
    
    est1(i) = pout.est1;
    est2(i) = pout.est2;
    
    % reconstruct F back from R,t
    T_rec = pout.est1.T;
    E_rec = skew(T_rec(1:3,4))*T_rec(1:3,1:3);
    F_rec = inv(params.K')*E_rec*inv(params.K);
    
    T_final = pout.est1.T_final;
    E_final = skew(T_final(1:3,4))*T_final(1:3,1:3);
    F_final = inv(params.K')*E_final*inv(params.K);
    
    pose_error = T_final\T_ss;
    r_err(1,i) = rot_error(pose_error);
    t_err(1,i) = trans_error(pose_error);
    
    % global optimization
    if i>3
        % converts feature indices into coordinates
        [tracks_x, tracks_y] = tracks_coords(info, M, 3, i);
        % figure;
        % imshow(i1,[]);
        % hold on;
        % plot(params.tracksx, params.tracksy);
        % compute cross-ratios        
        e = h2e(null(pout.est1.F));
        [ratio, sigma] = cross_ratio(tracks_x, tracks_y, e);
        stats(i).ratio = ratio;
        stats(i).sigma = sigma;

        x1 = [e2h(est2(i-1).x1(:,est2(i-1).inliers));...
              e2h(est2(i-1).x2(:,est2(i-1).inliers))];
        x2 = [e2h(est2(i).x1(:, est2(i).inliers));...
              e2h(est2(i).x2(:, est2(i).inliers))];

        x3 = [e2h(est1(i-1).x1); e2h(est1(i-1).x2)];
        x4 = [e2h(est1(i).x1); e2h(est1(i).x2)];
        
        % w is the weight of the cross ratio term in the optimization
        % objective
        w = linspace(.1, 5, 3);
        w = 1;
        for j = 1:length(w)
            t0  = [param.base 0 0]';
            x   = {x1, x2};
            T   = {inv(est2(i-1).T_final), inv(est2(i).T_final)};
            
            fun = @(c) objective1(w(j), param.K, t0, T, x, ratio, sigma, c);
            [c, ~, exitflag] = lsqnonlin(fun, [1, 1]);

            T = inv(est1(i).T_final);
            t = T(1:3, 4);
            t = c(2)*t;
            T(1:3, 4) = t;
            field = ['w1_',num2str(j)];
            if ~isfield(stats, field)
                stats(i).(field) = struct('c',[],...
                                       'c_prv',[],...
                                       'w',[],...
                                       'exitflag',[],...
                                       'T',[]);
                results{end+1} = field;
            end
            
            % update the results
            stats(i).(field).c = c;
            stats(i).(field).w = w(j);
            stats(i).(field).exitflag = exitflag;
            stats(i).(field).T = inv(T);
            
            continue;
            % operating points
            T  = {inv(est1(i-1).T_final), inv(est1(i).T_final)};            
            T1 = T{1}; R1 = T1(1:3,1:3); t1 = T1(1:3,4)';
            T2 = T{2}; R2 = T2(1:3,1:3); t2 = T2(1:3,4)';
            c0 = [0 0 0 t1 0 0 0 t2];
            h0(1,:) = quaternion.rotationmatrix(R1).e;
            h0(2,:) = quaternion.rotationmatrix(R2).e;
            x   = {x1, x2, x3, x4};
            fun = @(c) objective2(w(j), t0, x, ratio, sigma, h0, c);
            [c, ~, exitflag] = lsqnonlin(fun, c0);

            R = quaternion(param2quaternion(c(7:9)', h0(2,:)')).RotationMatrix;
            t = c(10:12)';
            T = [R t; 0 0 0 1];
            
            field = ['w2_',num2str(j)];
            if ~isfield(stats, field)
                stats(i).(field) = struct('c',[],...
                                       'c_prv',[],...
                                       'w',[],...
                                       'exitflag',[],...
                                       'T',[]);
                results{end+1} = field;
            end
            
            % update the results
            stats(i).(field).c = c;
            stats(i).(field).w = w(j);
            stats(i).(field).exitflag = exitflag;
            stats(i).(field).T = inv(T);
            
        end
    end
    
    % no optimisation result
    field = 'no_opt';
    if ~isfield(stats, field)
        stats(i).(field) = struct('T',[]);
    end
    stats(i).(field).T = est1(i).T_final;
    
    % DEBUG
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
    
    pin.F = {F_final, pout.est1.F, F_rec, F_gt, F_ss};
    pin.E = {E_final, pout.est1.E, E_rec, E_gt, E_ss};
    pin.T = {T_final, pout.est1.T, T_rec, T_gt, T_ss};
    
    pin.name = {'final T', 'epipolar estimation', 'after reconstruction', 'gt', 'ss'};
    
    %test_est_F(pin);
    
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
    
    %     test_est_F(pin);
end

% save results
for i=1:length(results)
    field = results{i};
    poses = nan(4, 4, num_frames);
    poses(:,:,1) = inv([eye(3) zeros(3,1); 0 0 0 1]);
    for j=2:num_frames
        if isempty(stats(j).(field))
            T = stats(j).no_opt.T;
        else
            T = stats(j).(field).T;
        end
        
        poses(:, :, j) = poses(:,:,j-1)/T;
    end
    
    dir = fullfile(RESULT_DIR, field);
    if exist(dir,'dir')
        command = ['rm -fr ', dir];
        system(command);
    end
    command = ['mkdir -p ', fullfile(dir, 'data')];
    system(command);
    savePoses(fullfile(dir, 'data', [sequence, '.txt']), poses);
end

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

function [T1,T2] = xr_obj(K,R,t,T1,T2,xr,x1,x2,x3,x4)

E = skew(t)*T;
F = inv(K')*E*inv(K);

dist_sampson(F,x1,x2);
end

function mt = complete_circle(mt, m12)

mt(4, :) = nan;
for j = 1:length(mt)
    ind = find(m12(1, :) == mt(1, j));
    if ~isempty(ind)
        mt(4, j) = m12(2, ind);
    end
end
end

function X = triangulate_points(c1p_ss, c2p_ss, num_pts, param)
X = nan(4, num_pts);
for j = 1:num_pts
    X(:, j) = vgg_X_from_xP_nonlin([c1p_ss(:, j) c2p_ss(:,j)], {param.P1, param.P2});
end
X = h2e(X);
end

