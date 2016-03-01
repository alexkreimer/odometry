function run_sequence(sequence, varargin)
close all
dbstop if error;

DATA_ROOT  = '/home/kreimer/KITTI/';
KITTI_HOME = fullfile(DATA_ROOT, 'dataset');

image_dir  = fullfile(KITTI_HOME, 'sequences', sequence);
poses_file = fullfile(KITTI_HOME, 'poses', [sequence, '.txt']);

D = dir(fullfile(image_dir, 'image_0','*.png'));
default_last = length(D(not([D.isdir])))-1;

p = inputParser;
p.addOptional('first',2,@isnumeric);
p.addOptional('last',default_last,@isnumeric);
p.addOptional('depth_thr',100,@isnumeric);
p.addOptional('inlier_thr',1,@isnumeric);
p.addOptional('ransac_iter',2,@isnumeric);
p.addOptional('sha','');
p.addOptional('run_ss',true,@isboolean);

p.KeepUnmatched = true;
parse(p,varargin{:});

RESULT_DIR = fullfile(DATA_ROOT, sprintf('results_%d_%d_%d',p.Results.depth_thr,p.Results.inlier_thr,p.Results.ransac_iter));

% setup camera parameters (KITTI)
[P0, P1] = util.kitti_read_calib(image_dir);
poses_gt = util.read_poses(poses_file);
param    = kitti_params(P0, P1);
K = param.K;

gt = process_gt(poses_gt);

poses1 = nan(4, 4, p.Results.last);
poses1(:,:,1) = inv([eye(3) zeros(3,1); 0 0 0 1]);
poses2 = poses1;

tracks = containers.Map('KeyType','uint64','ValueType','any');
update_tracks(fullfile('/home/kreimer/KITTI/tracks',sequence, sprintf('%06d.txt',p.Results.first-2)), tracks);

for i = p.Results.first:p.Results.last
    fprintf('processing frame %d of %d\n', i, p.Results.last);
    
    [i1, i2] = read_images(image_dir, i);
    [i1p,i2p]= read_images(image_dir, i-1);
    
    update_tracks(fullfile('/home/kreimer/KITTI/tracks',sequence,sprintf('%06d.txt',i-1)), tracks);
    
    coords2 = get_track_coords(tracks,2);
    
    x1p = permute(coords2(1:2,2,:), [1 3 2]);
    x2p = permute(coords2(3:4,2,:), [1 3 2]);
    [Xp,visiblep] = util.triangulate_chieral(x1p,x2p,param.P1,param.P2);
    
    x1 = permute(coords2(1:2,1,:), [1 3 2]);
    x2 = permute(coords2(3:4,1,:), [1 3 2]);
    [X,visible] = util.triangulate_chieral(x1,x2,param.P1,param.P2);

    directory = fullfile(DATA_ROOT,'data',sequence);
    filename = fullfile(directory, [sprintf('%06d',i),'.mat']);
    if ~exist(directory,'dir')
        command = ['mkdir -p ', fullfile(directory)];
        system(command);
    end    
    save(filename,'x1','x2','x1p','x2p','X','Xp','visible','visiblep');
    % use visibility constraint to filter wrong matches
    visible = visible&visiblep;
    x1  = x1(:,visible);
    x2  = x2(:,visible);
    x1p = x1p(:,visible);
    x2p = x2p(:,visible);
    X   = X(:,visible);
    Xp  = Xp(:,visible);
    
    if p.Results.run_ss
        disp('original reprojection minimization')
        tic;[a_ss, inliers,residual] = estimation.ransac_minimize_reproj(Xp, [x1; x2], param);toc;
        stats(i).ss.T                = util.tr2mat(a_ss);
        stats(i).ss.inliers          = {inliers};
        stats(i).ss.residual         = {residual};
    end
    
    %     disp('new reprojection minimization')
    %     tic; pose = estimation.ransac_min_reproj(K,param.base,Xp,x1,x2); toc;
    %     stats(i).ss1.T = pose;
    
    % These methods estimate R by decomposing F and then estimate t
    % separately
%     disp('F estimation/decomposition');
%     tic;
%     [TF1,F1,inliers1,residual1,success1] = estimation.rel_motion_F(K,x1p,x1);
%     [TF2,F2,inliers2,residual2,success2] = estimation.rel_motion_F(K,x2p,x1);
%     TF = estimation.stereo_motion_triangulate(TF1,TF2,[param.base 0 0]');
%     toc;
%     stats(i).F.T        = TF;
%     stats(i).F.inliers  = {inliers1,inliers2};
%     stats(i).F.residual = {residual1,residual2};
%     stats(i).F.success  = {success1,success2};
    
    disp('H_inf estimation');
    tic;
    depth = X(3,:);
    if sum(depth>50)>=10
        T1 = stats(i).ss.T;
        R = T1(1:3,1:3);
        t1 = T1(1:3,4);
        t2 = t1+R'*[param.base 0 0]';
        F1 = K'\util.skew(t1)*R/K;
        F2 = K'\util.skew(t2)*R/K;
        [TH1,inliers1,residual1] = estimation.rel_motion_H(K,[x1p x2p],[x1 x2],...
            [X(3,:) X(3,:)],param.base,'F',F1,'absRotInit',true,...
            'depth_thr',p.Results.depth_thr,'inlier_thr',p.Results.inlier_thr,'ransac_iter',p.Results.ransac_iter);
%         [TH2,inliers2,residual2] = estimation.rel_motion_H(K,x2p,x1,X(3,:),param.base,'F',F2,'absRotInit',true,...
%             'depth_thr',p.Results.depth_thr,'inlier_thr',p.Results.inlier_thr,'ransac_iter',p.Results.ransac_iter);
%         TH = estimation.stereo_motion_triangulate(TH1,TH2,[param.base 0 0]');
%         stats(i).Hg.sucess   = true;
%         stats(i).Hg.T        = TH;
%         stats(i).Hg.inliers  = {inliers1,inliers2};
%         stats(i).Hg.residual = {residual1,residual2};
        stats(i).Hg.success = true;
        stats(i).Hg.T = TH1;
        stats(i).Hg.inliers = inliers1;
        stats(i).Hg.residual = residual1;
    else
        if p.results.run_ss
            stats(i).Hg.T      = stats(i).ss.T;
            stats(i).Hg.sucess = false;
        else
            [a_ss,~,~] = estimation.ransac_minimize_reproj(Xp, [x1; x2], param);
            stats(i).Hg.T = util.tr2mat(a_ss);
            stats(i).Hg.sucess = false;
        end
    toc;
    
    stats(i).Hs.T        = stats(i).Hg.T;
    stats(i).Hs.T(1:3,4) = stats(i).ss.T(1:3,4);
    
    disp('One point translation estimation');
    R = stats(i).Hg.T(1:3,1:3);
    tic;[t,~,inliers] = estimation.ransac_minimize_reproj1(Xp,R,x1,x2,param);toc;
    stats(i).HX.T = [R t; 0 0 0 1];    
    stats(i).HX.inliers = inliers;
    
    st = stats(i);
    save(filename,'st','-append');
    
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
    
    % convert all poses to be relative to the world origin and save
    fields = fieldnames(stats);
    for ii=1:length(fields)
        field = fields{ii};
        poses = nan(4, 4, p.Results.last);
        poses(:,:,1) = [eye(3) zeros(3,1); 0 0 0 1];
        for j=2:i
            TH = stats(j).(field).T;
            poses(:,:,j) = poses(:,:,j-1)/TH;
        end
        directory = fullfile(RESULT_DIR, field);
        filename = fullfile(directory, 'data', [sequence, '.txt']);
        if ~exist(directory,'dir')
            command = ['mkdir -p ', fullfile(directory, 'data')];
            system(command);
        end
        util.savePoses(filename, poses);
    end
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

