function run_sequence(sequence, begin, num_frames)

close all
dbstop if error;

DATA_ROOT  = '/media/kreimer/my_drive/KITTI/';
KITTI_HOME = fullfile(DATA_ROOT, 'dataset');
RESULT_DIR = fullfile(DATA_ROOT, 'results');
DBG_DIR    = fullfile(DATA_ROOT, 'debug');

image_dir  = fullfile(KITTI_HOME, 'sequences', sequence);
poses_file = fullfile(KITTI_HOME, 'poses', [sequence, '.txt']);

if nargin<2
    begin = 2;
    num_frames = 10;
end

% setup camera parameters (KITTI)
[P0, P1] = util.kitti_read_calib(image_dir);
poses_gt = util.read_poses(poses_file);
param    = kitti_params(P0, P1);
K = param.K;

gt = process_gt(poses_gt);

poses1 = nan(4, 4, num_frames);
poses1(:,:,1) = inv([eye(3) zeros(3,1); 0 0 0 1]);
poses2 = poses1;

tracks = containers.Map('KeyType','uint64','ValueType','any');
update_tracks(fullfile('/media/kreimer/my_drive/KITTI/dataset/tracks',sequence, sprintf('%06d.txt',begin-2)), tracks);

for i = begin:num_frames
    fprintf('processing frame %d\n', i);
    
    [i1, i2] = read_images(image_dir, i);
    [i1p,i2p]= read_images(image_dir, i-1);

<<<<<<< HEAD
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

    % 3d based reprojection error minimization
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
    
    pout = stereo_motion_F(params);
=======
    update_tracks(fullfile('/media/kreimer/my_drive/KITTI/dataset/tracks',sequence,sprintf('%06d.txt',i-1)), tracks);
    
    coords2 = get_track_coords(tracks,2);
>>>>>>> 496c88ea7e05a7099c6cabaa46b92d2f7d05ddcf
    
    x1p = permute(coords2(1:2,2,:), [1 3 2]);
    x2p = permute(coords2(3:4,2,:), [1 3 2]);
    [Xp,visible] = triangulate_chieral(x1p,x2p,param.P1,param.P2);
    
    x1 = permute(coords2(1:2,1,:), [1 3 2]);
    x2 = permute(coords2(3:4,1,:), [1 3 2]);
    
    [X,~] = triangulate_chieral(x1,x2,param.P1,param.P2);
    
    Xp = Xp(:,visible);
    X  = X(:,visible);
    x1 = x1(:,visible);
    x2 = x2(:,visible);
    x1p= x1p(:,visible);
    x2p= x2p(:,visible);
    
    [a_ss, ~, ~, ~, ~] = estimation.ransac_minimize_reproj(Xp, [x1; x2], param);
    stats(i).ss.T = util.tr2mat(a_ss);

    %util.plot_triangles(i1,i2,i1p,i2p,x1(:,~visible),x2(:,~visible),x1p(:,~visible),x2p(:,~visible));
    %util.plot_triangles(i1,i2,i1p,i2p,x1(:,visible),x2(:,visible),x1p(:,visible),x2p(:,visible));
   
    [TF1,~] = estimation.rel_motion_F(K,x1p,x1);
    [TF2,inliers] = estimation.rel_motion_F(K,x2p,x1);
    TF = estimation.stereo_motion_triangulate(TF1,TF2,[param.base 0 0]');
    
    stats(i).no_opt_F.T = TF1;
    
    % keep this to use cross ratio later
    stats(i).no_opt_F.T2 = TF2;
    stats(i).no_opt_F.x1 = x2p(:,inliers);
    stats(i).no_opt_F.x2 = x1(:,inliers);

    TH1 = estimation.rel_motion_H(K,x1p,x1,X(3,:),param.base);
    TH2 = estimation.rel_motion_H(K,x2p,x1,X(3,:),param.base);
    TH = estimation.stereo_motion_triangulate(TH1,TH2,[param.base 0 0]');
    stats(i).no_opt_H.T = TH;
    
    R = estimation.H_inf_nonlin(K,x1,x2,X(3,:),param.base,150,.01);
    TX = estimation.trans_X(K,R,param.base,Xp,x1,x2);
    stats(i).tx.T = TX;
    
    % refinement
    if 0
        if i>3
            coords3 = get_track_coords(tracks,3);
            e = null(F);
            [ratio, sigma] = estimation.cross_ratio(coords3, e);
            
            tr1.T = stats(i).no_opt_F.T2;
            tr1.x1= stats(i).no_opt_F.x1;
            tr1.x2= stats(i).no_opt_F.x2;
            
            tr2.T = stats(i-1).no_opt_F.T2;
            tr2.x1= stats(i-1).no_opt_F.x1;
            tr2.x2= stats(i-1).no_opt_F.x2;
            
            fun = @(c) estimation.objective1(1,K, param.base,...
                tr1, tr2, ratio, sigma, c);
            [c, ~, exitflag] = lsqnonlin(fun, [1 1]);
            
            T = inv(stats(i).no_opt_F.T);
            T(1:3,4) = c(2)*T(1:3,4);
            field = ['objective1_w_',strrep(num2str(w),'.','')];
            %update the results
            stats(i).(field).T = inv(T);
            %
            %             continue;
            %             %operating points
            %             T  = {inv(est1(i-1).T_final), inv(est1(i).T_final)};
            %             T1 = T{1}; R1 = T1(1:3,1:3); t1 = T1(1:3,4)';
            %             T2 = T{2}; R2 = T2(1:3,1:3); t2 = T2(1:3,4)';
            %             c0 = [0 0 0 t1 0 0 0 t2];
            %             h0(1,:) = quaternion.rotationmatrix(R1).e;
            %             h0(2,:) = quaternion.rotationmatrix(R2).e;
            %             x   = {x1, x2, x3, x4};
            %             fun = @(c) objective2(w(j), t0, x, ratio, sigma, h0, c);
            %             [c, ~, exitflag] = lsqnonlin(fun, c0);
            %
            %             R = quaternion(param2quaternion(c(7:9)', h0(2,:)')).RotationMatrix;
            %             t = c(10:12)';
            %             T = [R t; 0 0 0 1];
            %
            %             field = ['w2_',num2str(j)];
            %             if ~isfield(stats, field)
            %                 stats(i).(field) = struct('c',[],...
            %                     'c_prv',[],...
            %                     'w',[],...
            %                     'exitflag',[],...
            %                     'T',[]);
            %                 results{end+1} = field;
            %             end
            %
            %             update the results
            %             stats(i).(field).c = c;
            %             stats(i).(field).w = w(j);
            %             stats(i).(field).exitflag = exitflag;
            %             stats(i).(field).T = inv(T);
            %         end
        else
            for j=1:1
                field = ['objective1_w_',strrep(num2str(w),'.','')];
                stats(i).(field).T = stats(i).no_opt_F.T;
            end
        end
    end
    
end

% convert all poses to be relative to the world origin and save
fields = fieldnames(stats);
for i=1:length(fields)
    field = fields{i};
    poses = nan(4, 4, num_frames);
    poses(:,:,1) = [eye(3) zeros(3,1); 0 0 0 1];
    for j=2:num_frames
        dist(j) = norm(stats(j).(field).T(1:3,4));
    end
    dist
    for j=2:num_frames
        TH = stats(j).(field).T;
        poses(:,:,j) = poses(:,:,j-1)/TH;
    end
    for j=2:num_frames
        dist_after(j) = norm(poses(1:3,4,j-1)-poses(1:3,4,j));
    end
    dist_after
    dir = fullfile(RESULT_DIR, field);
    if exist(dir,'dir')
        command = ['rm -fr ', dir];
        system(command);
    end
    command = ['mkdir -p ', fullfile(dir, 'data')];
    system(command);
    util.savePoses(fullfile(dir, 'data', [sequence, '.txt']), poses);
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

function [X,visible] = triangulate_chieral(x1,x2,P1,P2,i1,i2)
num_pts = length(x1);
X = nan(4,num_pts);
visible = false(1,num_pts);
detM1 = det(P1(1:3,1:3));
detM2 = det(P2(1:3,1:3));
for j = 1:num_pts
    if nargin > 4
    figure;
    imshow([i1;i2]); hold on; plot([x1(1,j);x2(1,j)],[x1(2,j);x2(2,j)+size(i1,1)]);
    end
    
    X(:, j) = vgg.vgg_X_from_xP_nonlin([x1(:, j) x2(:,j)],{P1,P2});
    
    p1 = P1*X(:,j);
    p2 = P2*X(:,j);
    visible(j) = X(4,j)*p1(3)*detM1 > 0 && X(4,j)*p2(3)*detM2 > 0;
    x1(:,j)'
    x2(:,j)'
    visible(j)
end
X = util.h2e(X);
end

function X = triangulate_points(c1p_ss, c2p_ss, param)
num_pts = length(c1p_ss);
X = nan(4, num_pts);
for j = 1:num_pts
    X(:, j) = vgg.vgg_X_from_xP_nonlin([c1p_ss(:, j) c2p_ss(:,j)], {param.P1, param.P2});
end
X = util.h2e(X);
end

function update_tracks(file_name, tracks)
fd = fopen(file_name,'r');
A = fscanf(fd, '%d %d %d %d %d',[5 Inf]);
fclose(fd);

k1 = cell2mat(keys(tracks));
k2 = A(1,:);

delta = setdiff(k1,k2);
remove(tracks, num2cell(delta));

for i=1:size(A,2)
    tid = A(1,i);
    tval= A(2:end,i);
    if isKey(tracks,tid)
        val = tracks(tid);
        tracks(tid) = [tval val];
    else
        tracks(tid) = tval;
    end
end

end

function coords =  get_track_coords(tracks, len)

k = keys(tracks);
t = 1;
coords = nan(4,len,0);

for i=1:length(k)
    val = tracks(k{i});
    if size(val,2)>=len
        coords(:,1:len,t) = val(:,1:len);
        t = t+1;
    end
end

end

