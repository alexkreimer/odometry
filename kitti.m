function kitti()
close all;
dbstop if error;

KITTI_HOME = '/home/kreimer/tmp/KITTI/dataset';

image_dir = fullfile(KITTI_HOME, 'sequences', '01');
poses_file = fullfile(KITTI_HOME, 'poses','01.txt');

first = 1;
last = 5;

% setup camera parameters (KITTI)
[P0, P1] = kitti_read_calib(image_dir);
gt_poses = kitti_read_poses(poses_file);

% baseline, focal, principal point
param.base = -P1(1,4)/P1(1,1);
param.calib.f = P0(1,1);
param.calib.cu = P0(1,3);
param.calib.cv = P0(2,3);
% this commands minimization procedure to first solve 3d-3d rigid motion
% and use it as an initial point for Gauss-Newton reprojection minimization
% procedure
param.init = true;
% fundamental
param.F = vgg_F_from_P(P0,P1);
% number of corners to extract per image
param.corner_num = 1000;
% descriptors are of size (2win_sz+1)x(2win_sz+1)
param.win_sz = 5;
% minimal acceptable disparity
param.min_d = 2;
% BA window
param.ba_w = 3;

% observation dimension
param.obd = 4;
% camera parameter vector dimension
param.ad = 6;
% structure parameter vector dimension
param.bd = 3;

param.lm_max_iter = 1000;

% init KLT trackers, 1 for the left cam, 2 for the right
klt1 = vision.PointTracker;
klt2 = vision.PointTracker;

a = repmat({zeros(6,1)},3,1);

% main loop, iteration per image pair
for i=first:last
    % read images
    [i1,i2] = read_kitti_images(image_dir, i);

    if i==first,
        ft1 = feat_get(i1, param);
        ft2 = feat_get(i2, param);
    else
        [ft1, valid1] = feat_track(i1,klt1,param);
        [ft2, valid2] = feat_track(i2,klt2,param);
    end

    % match features across the stereo pair, using epipolar constraint
    match = feat_match(i1,i2,ft1,ft2,'debug',false,'F',param.F,'do_sampson',true);
    
    %ft1 = ft1(match(1,:)); ft2 = ft2(match(2,:));
    len = length(ft1);
    
    if i==first
        tracks = track_init(len);
    else
        %ft1_p = ft1_p(match(1,:)); ft2_p = ft2_p(match(2,:));
    end

    
    % create/update tracks
    for j = 1:len
        % circle test
        % ind = find(match(1,:),j); if match(2,ind) ~= ind, continue, end;
        x1 = ft1(match(1,j)).pt;
        x2 = ft2(match(2,j)).pt;
        [X, d, J] = triangulate1(x1,x2,param);
        tracks(j) = track_update(tracks(j),x1,x2,X,d,J,i,param);
    end
    
    % debug
    if i>1
        do_debug(i1,i2,i1_p,i2_p,ft1,ft2,ft1_p,ft2_p,tracks,param);
    end
    
    % delete dead tracks
    ft1 = ft1([tracks.live]);
    ft2 = ft2([tracks.live]);
    if i>1
        ft1_p = ft1_p([tracks.live]);
        ft2_p = ft2_p([tracks.live]);
    end
    tracks([tracks.live]==false) = [];
    [x,X] = track_getx(tracks);
    if i==1
        % init the KLT trackers with matched features
        initialize(klt1, x(1:2,:)',i1);
        initialize(klt2, x(3:4,:)',i2);
        a{1} = zeros(6,1);
    else
        setPoints(klt1,x(1:2,:)');
        setPoints(klt2,x(3:4,:)');
    end
    
    % solve incremental rigid motion
    [a{i}, ~] = ransac_minimize_reproj(X,x,param);

    %[x,~] = track2cell(tracks);
    %if i>= size(a,1)
    %    sigma = diag(ones(1,4*numel(x)));
    %    [aba,bba] = simple_ba(x,sigma,a,b,param);
    %end
    i1_p = i1; i2_p = i2;
end
end


function [i1,i2] = read_kitti_images(seq_home, idx)

file_names = {fullfile(seq_home, 'image_0', sprintf('%06d.png',idx))
    fullfile(seq_home, 'image_1', sprintf('%06d.png',idx))};

i1 = imread(file_names{1});
i2 = imread(file_names{2});

end

function poses = kitti_read_poses(poses_file)
fid = fopen(poses_file, 'r');
i=1;
tline = fgetl(fid);
while ischar(tline)
   poses(:,:,i) = reshape(sscanf(tline, '%f %f %f %f %f %f %f %f %f %f %f %f'),[4 3])';
   tline = fgetl(fid);
   i=i+1;
end
fclose(fid);
end

function [P0, P1] = kitti_read_calib(seq_home)

calib_file = fullfile(seq_home, 'calib.txt');

fd = fopen(calib_file, 'r');
p0 = fgetl(fd);
p1 = fgetl(fd);
P0 = reshape(sscanf(p0, 'P0: %f %f %f %f %f %f %f %f %f %f %f %f'),[4 3])';
P1 = reshape(sscanf(p1, 'P1: %f %f %f %f %f %f %f %f %f %f %f %f'),[4 3])';
fclose(fd);
end

function [features, valid] = corners2features(i1,c,win_sz)
% corners may become invalid because they are too close to image border
pts_num = size(c,1);
features = struct('pt',cell(pts_num,1),'d',nan);
valid = true(pts_num,1);
for j=1:pts_num
    x1 = c(j,1)-win_sz;
    x2 = c(j,1)+win_sz;
    y1 = c(j,2)-win_sz;
    y2 = c(j,2)+win_sz;
    if x1<1 || y1<1 || x2>size(i1,2) || y2>size(i1,1)
        valid(j) = false;
        continue;
    else
        patch = i1(y1:y2,x1:x2);
        features(j).pt = c(j,:)';
        features(j).d = single(patch(:));
    end
end

end

function show4(i1,i2,i3,i4,x1,x2)
i = [[i1,i2];[i3,i4]];

imshow(i); hold on;
plot([x1(1) x1(3)+size(i1,2) x2(3)+size(i1,2) x2(1) x1(1)], [x1(2) x1(4) x2(4)+size(i1,1) x2(2)+size(i1,1) x1(2)]);
hold off;
end


function [x,X] = track2cell(track)
X = repmat({nan(3,1)},length(track),1);
for i=1:length(track)
    for j=1:length(track(i).ft)
        x{i,j} = track(i).ft(j).pts;
    end
    X{i} = track(i).ft(j).X;
end
end

function frames = feat_get(i1, param)
corners = corner(i1,'harris',param.corner_num);
[frames,valid] = corners2features(i1,round(corners),param.win_sz);
frames = frames(:,valid);
end


function [frames, validity] = feat_track(i1,klt,param)
[corners, validity] = step(klt,i1);
[frames,validity1] = corners2features(i1,round(corners),param.win_sz);
validity = validity & validity1;
end

function [X,d,J] = triangulate1(x1,x2,param)
[X, d, J] = triangulate_naive(x1,x2,param.base,param.calib.f,...
    param.calib.cu,param.calib.cv);
end

function track = track_update(track,x1,x2,X,d,J,i,param)
track.len = track.len+1;
track.ft(track.len).x = [x1;x2];
track.ft(track.len).X = X;
track.ft(track.len).J = J;
track.ft(track.len).d = d;
track.live = d>param.min_d;
track.updated_at = i;
end

function tracks = track_init(N)
s = struct('ft',struct('x',[],'X',nan,'J',nan,'d',nan),'live',nan,...
    'updated_at',nan,'len',0);
tracks = repmat(s,[N 1]);
end

function do_debug(i1,i2,i1_p,i2_p,ft1,ft2,ft1_p,ft2_p,tracks,param)

x(1:2,:) = [ft1.pt];
x(3:4,:) = [ft2.pt];
x_p(1:2,:) = [ft1_p.pt];
x_p(3:4,:) = [ft2_p.pt];

figure; showMatchedFeatures(i1,i1_p,x(1:2,:)',x_p(1:2,:)','blend'); title('left/prev');
figure; showMatchedFeatures(i2,i2_p,x(3:4,:)',x_p(3:4,:)','blend'); title('right/prev');

x = nan(4,length(tracks));
for j=1:length(tracks)
    x(:,j) = tracks(j).ft(end).x; d(j) = abs(tracks(j).ft(end).d);
end
figure; showMatchedFeatures(i1,i2,x(1:2,:)',x(3:4,:)','blend');
hold on; title('disparity map');
for j=1:length(tracks),
    text(x(1,j),x(2,j),sprintf('%.2g',param.base*param.calib.f/d(j)),'FontSize',8)
end
d(d==0) = .1;
scatter(x(1,:),x(2,:),d,'b'); legend('dead','live');
end

function [x,X] = track_getx(tracks)
len = length(tracks);
x = nan(4,len);
X = nan(3,len);
for j=1:len
    x(:,j) = tracks(j).ft(end).x;
    if nargout == 2
        X(:,j) = tracks(j).ft(1).X;
    end
end
end

