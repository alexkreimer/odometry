function run_sequence(sequence, varargin)
close all
dbstop if error;

DATA_ROOT       = '/home/kreimer/KITTI/';
%DATA_ROOT = '/media/kreimer/my_drive/odometry_data/';
KITTI_HOME      = fullfile(DATA_ROOT, 'dataset');
% minimal number of distant points; if there are less, run stereo scan
NUM_DIST_THRESH = 10;

% input images dir
image_dir       = fullfile(KITTI_HOME, 'sequences', sequence);

% distance threshold in meters
dist_thresh     = 50;

D = dir(fullfile(image_dir, 'image_0','*.png'));
default_last = length(D(not([D.isdir])))-1;

p = inputParser;
p.addOptional('first',2,@isnumeric);
p.addOptional('last',default_last,@isnumeric);
p.addOptional('depth_thr',100,@isnumeric);
p.addOptional('inlier_thr',1,@isnumeric);
p.addOptional('ransac_iter',2,@isnumeric);
p.addOptional('sha','');
p.addOptional('mono_left',0,@isnumeric);
p.addOptional('mono_right',0,@isnumeric);
p.addOptional('stereo',0,@isnumeric);
p.addOptional('Fp',0,@isnumeric);
p.addOptional('monomono',0,@isnumeric);

p.KeepUnmatched = true;
parse(p,varargin{:});

% output directories
RESULT_DIR      = fullfile(DATA_ROOT, 'info', sprintf('results_%s', p.Results.sha));
RESDATA_DIR     = fullfile(DATA_ROOT, 'info', sprintf(   'data_%s', p.Results.sha),sequence);

%if exist(RESDATA_DIR,'dir')
%    system(['rm -fr ', RESDATA_DIR]);
%end
system(['mkdir -p ', fullfile(RESDATA_DIR)]);

%if exist(RESULT_DIR, 'dir')
%    system(['rm -fr ', RESULT_DIR]);       
%end
system(['mkdir -p ', fullfile(RESULT_DIR)]);

C = cellfun(@num2str, struct2cell(p.Results), 'UniformOutput', false);
C(:,2) = fieldnames(p.Results);
T = cell2table(C,'VariableNames',{'Value','Name'});
writetable(T,fullfile(RESDATA_DIR, 'params.dat'))

% setup camera parameters (KITTI)
[P0, P1] = util.kitti_read_calib(image_dir);
param    = kitti_params(P0, P1);
K = param.K;

tracks = containers.Map('KeyType','uint64','ValueType','any');
update_tracks(fullfile(DATA_ROOT, 'tracks', sequence, sprintf('%06d.txt',p.Results.first-2)), tracks);

first = true;
for i = p.Results.first:p.Results.last
    fprintf('processing frame %d of %d\n', i, p.Results.last);
    
    [i1, i2] = read_images(image_dir, i);
    [i1p,i2p]= read_images(image_dir, i-1);

    update_tracks(fullfile(DATA_ROOT,'tracks',sequence,sprintf('%06d.txt',i-1)), tracks);
    coords2 = get_track_coords(tracks,2);
    if ~first
        coords3 = get_track_coords(tracks, 3);

        % tracks of length 3
        if size(coords3, 3)
            tr3.x_ppv = permute(coords3(1:2,3,:), [1 3 2]);
            tr3.x_prv = permute(coords3(1:2,2,:), [1 3 2]);
            tr3.x_cur = permute(coords3(1:2,1,:), [1 3 2]);
        else
            tr3.x_ppv = [];
        end
    end

    x1p = permute(coords2(1:2,2,:), [1 3 2]);
    x2p = permute(coords2(3:4,2,:), [1 3 2]);
    [Xp,visiblep] = util.triangulate_chieral(x1p,x2p,param.P1,param.P2);
    
    x1 = permute(coords2(1:2,1,:), [1 3 2]);
    x2 = permute(coords2(3:4,1,:), [1 3 2]);
    [X,visible] = util.triangulate_chieral(x1,x2,param.P1,param.P2);
    
    % save some stats
    filename = fullfile(RESDATA_DIR, 'stats.mat');
    save(filename,'x1','x2','x1p','x2p','X','Xp','visible','visiblep','-v7.3');
    
    % use visibility constraint to filter wrong matches
    visible = visible&visiblep;
    x1  = x1(:,visible);
    x2  = x2(:,visible);
    x1p = x1p(:,visible);
    x2p = x2p(:,visible);
    X   = X(:,visible);
    Xp  = Xp(:,visible);

    % in mono setting we can not rely on triangulation to filter outliers,
    % so we keep the original matches
    tr2.x1_prv = x1p;
    tr2.x2_prv = x2p;
    
    tr2.x1_cur = x1;
    tr2.x2_cur = x2;
    
    disp('algorithm: reprojection minimization')
    tic;[a_ss, inliers,residual] = estimation.ransac_minimize_reproj(Xp, [x1; x2], param);toc;
    stats(i).ss.T                = util.tr2mat(a_ss);
    stats(i).ss.inliers          = {inliers};
    stats(i).ss.residual         = {residual};
    
    disp('algorithm: IO');
    if p.Results.monomono
        if first
            x_prv = tr2.x1_prv;
            x_cur = tr2.x1_cur;
            [mask, num_dist] = choose_distant_stereo(Xp, dist_thresh);
        else
            % no depth information available
            x_prv = tr2.x1_prv;
            x_cur = tr2.x1_cur;
            [mask, num_dist]  = choose_distant_mono(x_prv, x_cur, Hp, dist_thresh);
        end
    else
        tic;
        [mask, num_dist] = choose_distant_stereo(Xp, dist_thresh);
        if p.Results.mono_left
            % mono emulation
            x_prv = x1p;
            x_cur = x1;
        elseif p.Results.mono_right
            x_prv = x2p;
            x_cur = x1;
        elseif p.Results.stereo
            x_prv = [x1p x2p];
            x_cur = [x1 x2];
        else
            error('you need to specify either mono_left or mono_right or stereo');
        end
    end
    
    if num_dist >= NUM_DIST_THRESH
        % SS computed fundamental
        F1 = trans2fund(stats(i).ss.T, K);

        % choose fundametal to use for the rotation estimation
        if p.Results.Fp
            if first
                F = F1;
            else
                F = Fp;
            end
        else
            F = F1;
        end

        % estimate rotation
        [Hp, R, ~, ~] = estimation.H_inf_nonlin(K,x_prv,x_cur,mask,'F',F,'absRotInit',true,...
                    'inlier_thr',p.Results.inlier_thr,'ransac_iter',p.Results.ransac_iter);

        t = estimation.trans_geom(K, Hp, x_prv, x_cur);
        t = t*norm(stats(i).ss.T(1:3,4));
        stats(i).HG.T = [R -R*t; 0 0 0 1];
        stats(i).HG.success = true;
        
        if p.Results.monomono
            if first
                t = stats(i).ss.T(1:3,4);
                %t = normc(t);
            else
                % previous motion estimates
                Tp = stats(i-1).HX.T;
                
                Rp = Tp(1:3,1:3);
                tp = Tp(1:3,4);
                
                % compose camera matrices
                P1 = K*[eye(3) zeros(3,1)];
                P2 = K*[Rp tp];
                
                % triangulate points
                [Xp, vis] = util.triangulate_chieral(tr3.x_ppv,tr3.x_prv,P1,P2);
                
                %imshow(i1p); hold on;
                %x_prj = util.h2e(P2*util.e2h(Xp(:,vis)));
                %plot([x_prj(1,:);tr3.x_prv(1,vis)], [x_prj(2,:);tr3.x_prv(2,vis)]);
                
                %showMatchedFeatures(i1,i1p,tr3.x_cur(:,vis(1:20))', tr3.x_prv(:,vis(1:20))', 'montage', 'PlotOptions', {'ro','g+','y-'});
                
                % minimize reprojection errors
                x_cur = tr3.x_cur(:,vis);
                [t, predict, inliers] = estimation.ransac_minimize_reproj1(Xp(:,vis), R, x_cur, param);
                T1 = [R t; 0 0 0 1];
                T2 = [Rp tp; 0 0 0 1];
                T = inv(T2)*T1;
                t = T(1:3,4);
                
                % plot predicted vs. observed
                imshow(i1);
                hold on; title('predicted vs. observed for HX');
                plot([x_cur(1,inliers); predict(1,inliers)], [x_cur(2,inliers); predict(2,inliers)], 'LineWidth', 2);
                
                % residuals histogram
                figure;
                dx = util.colnorm(predict(:,inliers) - x_cur(:, inliers));
                hist(dx);
            end
        else
            [t,~,inliers] = estimation.ransac_minimize_reproj1(Xp,R,[x1;x2],param);
        end
        
        stats(i).HX.T = [R t; 0 0 0 1];
        stats(i).HX.inliers = inliers;
        stats(i).HX.success = true;
    else
        % fall-back to HX
        stats(i).HX.T      = stats(i).ss.T;
        stats(i).HX.success= false;
        
        stats(i).HG.T      = stats(i).ss.T;
        stats(i).HG.success= false;        
    end
    
    t = stats(i).HX.T(1:3,4);
    R = stats(i).HX.T(3:3,3:3);
    Fp = K'\util.skew(t)*R/K;
    
    toc;
    save(filename,'stats','-append');

    % convert all poses to be relative to the world origin and save
    fields = fieldnames(stats);
    for ii=1:length(fields)
        field = fields{ii};
        poses = nan(4, 4, i);
        poses(:,:,1) = [eye(3) zeros(3,1); 0 0 0 1];
        for j=2:i
            pose = stats(j).(field).T;
            poses(:,:,j) = poses(:,:,j-1)/pose;
        end        
        directory = fullfile(RESULT_DIR, field);
        filename = fullfile(directory, 'data', [sequence, '.txt']);
        if ~exist(directory,'dir')
            command = ['mkdir -p ', fullfile(directory, 'data')];
            system(command);
        end
        util.savePoses(filename, poses);
    end
    
    first = false;
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

function F = trans2fund(T, K)
% compose the Fundamental matrix from the transformation parameters

T1 = T;
R  = T1(1:3,1:3);
t1 = T1(1:3,4);
F  = K'\util.skew(t1)*R/K;
end


function [mask, num_dist] = choose_distant_mono(x_prv, x_cur, Hp, thresh)
    Hx_prv = util.h2e(Hp*util.e2h(x_prv));
    
    delta = util.colnorm(Hx_prv-x_cur);
    mask = delta<thresh;
    
    if nargout > 1
        num_dist = sum(mask);
    end
end

function [mask, num_dist] = choose_distant_stereo(X, thresh)
    mask = X(3,:) > thresh;
    if nargout > 1
        num_dist = sum(mask);
    end
    
end
