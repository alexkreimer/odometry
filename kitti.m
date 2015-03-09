function kitti()
close all;
dbstop if error;

KITTI_HOME = '/home/kreimer/tmp/KITTI/dataset';


seq_home = fullfile(KITTI_HOME, 'sequences', '01');
first = 1;
last = 5;

% setup camera parameters (KITTI)
[P0, P1] = read_kitti_calib(seq_home);

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

% init KLT trackers, 1 for the left cam, 2 for the right
klt1 = vision.PointTracker;
klt2 = vision.PointTracker;

a = repmat({zeros(6,1)},3,1);

% main loop, iteration per image pair
for i=first:last
    % read images
    [i1,i2] = read_kitti_images(seq_home, i);
    
    % first pair is special, since we only init things
    if i==first,
        % extract harris corners and convert them to struct arrays of
        % features
        f1 = corners2features(i1,corner(i1,'harris',param.corner_num),param.win_sz);
        f2 = corners2features(i2,corner(i2,'harris',param.corner_num),param.win_sz);
        
        % match features across the stereo pair, using epipolar constraint
        match = match_features(i1,i2,f1,f2,'debug',false,'F',param.F,'do_sampson',true);
        match_num = size(match,2);
        
        % preallocate
        % track is the tracklets db
        track = struct('ft',cell(match_num,1));
        
        for j = 1:match_num
            % triangulate matches, returns 3d point, disparity and the Jacobian
            % of the triangulation func, which may be used to compute the
            % cov of this 3d point (using sandwitch cov formula)
            [X ,d, J] = triangulate_naive([f1(match(1,j)).pt],[f2(match(2,j)).pt],...
                param.base,param.calib.f,param.calib.cu,param.calib.cv);
            
            % store feature in the db
            track(j).ft(1).pts = [f1(match(1,j)).pt;f2(match(2,j)).pt];
            track(j).ft(1).X = X;
            track(j).ft(1).J = J;
            track(j).ft(1).d = d;
            track(j).begin = i;
            track(j).updated_at = i;
            track(j).live = abs(d)>param.min_d;
        end
        
        % delete dead tracks
        dead_tracks = track([track.live]==false);
        track([track.live]==false) = [];
        
        % debug
        xs = nan(2,length(dead_tracks)); dead_normJ = nan(length(dead_tracks),1);
        for j=1:length(dead_tracks), 
            xs(1,j) = dead_tracks(j).ft(1).pts(1); xs(2,j) = dead_tracks(j).ft(1).pts(2);
            dead_normJ(j) = norm(dead_tracks(j).ft(1).J*dead_tracks(j).ft(1).J');
        end
        figure; subplot(211); imshow(i1); hold on; plot(xs(1,:),xs(2,:),'xr');
        
        normJ = nan(length(track),1);
        for j=1:length(track), 
            xs(1,j) = track(j).ft(1).pts(1); xs(2,j) = track(j).ft(1).pts(2);
            normJ(j) = norm(track(j).ft(1).J*track(j).ft(1).J');
            d(j) = track(j).ft(1).d;
        end
        plot(xs(1,:),xs(2,:),'xg'); legend('dead','live');
        [normJ,idx] = sort(normJ); d = d(idx);
        subplot(212); plot(normJ,'b'); hold on; plot(abs(d),'g');

        % init the KLT trackers with matched features
        [x,~] = track2cell(track);
        pts = cell2mat(x(:,1)');
        initialize(klt1, pts(1:2,:)',i1);
        initialize(klt2, pts(3:4,:)',i2);
        
        % store current images (4debug)
        i1_prev = i1; i2_prev = i2;
        a{1} = zeros(6,1);
        continue; 
    end
    
    % track the corners into the current images and convert them to
    % features
    [points1, point_validity1] = step(klt1, i1);
    [points2, point_validity2] = step(klt2, i2);
    f1 = corners2features(i1,round(points1),param.win_sz);
    f2 = corners2features(i2,round(points2),param.win_sz);
    
    % match across the stereo pair
    match = match_features(i1,i2,f1,f2,'debug',false,'F',param.F,'do_sampson',true);
    
    % find the tracklets that survived current step
    for j=1:length(point_validity1)
        track(j).live = false;
        
        % KLT says we lost it...
        if ~point_validity1(j) || ~point_validity2(j), continue; end
        
        % circle test
        ind = find(match(1,:),j); if match(2,ind) ~= ind, continue, end;
        
        [X, d, J] = triangulate_naive(points1(j,:),points2(j,:),...
            param.base,param.calib.f,param.calib.cu,param.calib.cv);
        
        if d<param.min_d, continue; end  % TBD: maybe zero disparity is not a reason to kill this track?
        
        cur_track_len = length([track(j).ft]);
        track(j).ft(cur_track_len+1).pts = [points1(j,:)';points2(j,:)'];
        track(j).ft(cur_track_len+1).X = X;
        track(j).ft(cur_track_len+1).J = J;
        track(j).ft(cur_track_len+1).d = d;
        track(j).live = true;
        track(j).updated_at = i;
    end
    
    to_remove = false(length(track),1);
    for j=1:length(track)
        if track(j).live == false && (track(j).updated_at + param.ba_w < i)
            to_remove(j) = true;
        end
    end
    
    % remove the tracklets that got out of the BA window
    track(to_remove) = [];
    [x,b] = track2cell(track);
    pts = nan(4,length(track));
    k=1;
    for j=1:length(track)
        pts(:,j) = track(j).ft(end).pts;
        if track(j).live,
            X(:,k) = track(k).ft(1).X;
            observed(:,k) = track(k).ft(end).pts;
            k = k+1;
        end
    end
    
    % keep trackers in sync
    setPoints(klt1,pts(1:2,:)',[track.live]);
    setPoints(klt2,pts(3:4,:)',[track.live]);
    
    % solve rigid motion
    [a{i}, inliers] = ransac_minimize_reproj(X, observed, param);
    if i>= size(a,1)
        sigma = diag(ones(1,4*numel(x)));
        [aba,bba] = simple_ba(x,sigma,a,b,param);
    end
end

end

function [i1,i2] = read_kitti_images(seq_home, idx)

file_names = {fullfile(seq_home, 'image_0', sprintf('%06d.png',idx))
    fullfile(seq_home, 'image_1', sprintf('%06d.png',idx))};

i1 = imread(file_names{1});
i2 = imread(file_names{2});

end

function [P0, P1] = read_kitti_calib(seq_home)

calib_file = fullfile(seq_home, 'calib.txt');

fd = fopen(calib_file, 'r');
p0 = fgetl(fd);
p1 = fgetl(fd);
P0 = reshape(sscanf(p0, 'P0: %f %f %f %f %f %f %f %f %f %f %f %f'),[4 3])';
P1 = reshape(sscanf(p1, 'P1: %f %f %f %f %f %f %f %f %f %f %f %f'),[4 3])';
end

function ft = corners2features(i1,c,win_sz)
pts_num = size(c,1);
%ft = struct('pt',cell(pts_num,1));
j=1;
for i=1:pts_num
    x1 = c(i,1)-win_sz;
    x2 = c(i,1)+win_sz;
    y1 = c(i,2)-win_sz;
    y2 = c(i,2)+win_sz;
    if x1<1 || y1<1 || x2>size(i1,2) || y2>size(i1,1)
        continue;
    end
    patch = i1(y1:y2,x1:x2);
    ft(j).pt = c(i,:)';
    ft(j).d = single(patch(:));
    j = j+1;
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