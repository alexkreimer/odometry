function kitti_oo()
close all;
dbstop if error;

KITTI_HOME = '/home/kreimer/tmp/KITTI/dataset';
DBG_DIR = 'debug';

image_dir = fullfile(KITTI_HOME, 'sequences', '01');
poses_file = fullfile(KITTI_HOME, 'poses','01.txt');

first = 1;
last = 50;

% setup camera parameters (KITTI)
[P0, P1] = kitti_read_calib(image_dir);
poses = kitti_read_poses(poses_file);

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

param.feat_num = 500;                        % number of corners to extract per image
param.feat_method = 'MinimumEigenvalue';     % corner extraction method (not used)
param.patchr = 3;                            % descriptors are of size (2*patchr+1)**2
param.searchr = 100;                         % search window radius for feature tracking

param.min_d = 4;                             % minimal acceptable disparity
param.min_disp = 10;
param.ba_w = 3;                              % BA window

% reprojection error minimization error threshold
param.inlier_thresh = 10;
param.model_size = 3;
param.ransac_iter = 100;

% observation dimension
param.obd = 4;
% camera parameter vector dimension
param.ad = 6;
% structure parameter vector dimension
param.bd = 3;

param.lm_max_iter = 1000;


% init KLT trackers, 1 for the left cam, 2 for the right
klt1 = vision.PointTracker('BlockSize',[41 41],'MaxBidirectionalError',5,'MaxIterations',50);
klt2 = vision.PointTracker('BlockSize',[41 41],'MaxBidirectionalError',5,'MaxIterations',50);

a = repmat({zeros(6,1)},3,1);

% main loop, iteration per image pair
for i=first:last
    % read images
    [i1,i2] = read_kitti_images(image_dir, i);
    
    if i==first,
        features1 = features_init(i1,param);
        features2 = features_init(i2,param);
        tracklets1 = tracklets_init(features1,i);
        tracklets2 = tracklets_init(features2,i);
    else
        tic;features1 = tracklets_track(tracklets1,i1,param.searchr,param.patchr);toc;
        features2 = tracklets_track(tracklets2,i2,param.searchr,param.patchr);
        
        % track features independently in the left and right views
        %features1 = features_track(i1,klt1,param);
        %features2 = features_track(i2,klt2,param);

        % note that for these the indices of the features do not correspond
        % the indices of the tracklets (this is because invalids are
        % discarded)
        pts1 = features_pts(features1);
        pts2 = features_pts(features2);
        % update tracklets.  Tracks that are lost are not updated and thus
        % their valid() call will be false
        tracklets1 = tracklets_update(tracklets1,features1,features_augment(i1,pts1,param),i);
        tracklets2 = tracklets_update(tracklets2,features2,features_augment(i1,pts2,param),i);
        
        % plot left/right tracks
        tracklets_plot(tracklets1,i1,pi1,sprintf('frame %d: left vs prev',i));
        tracklets_plot(tracklets2,i2,pi2,sprintf('frame %d: right vs prev',i));
    end

    
    % match between current left/right pair
    match = horiz_match(tracklets1,tracklets2,param);
    
    % plot pair matches
    %tracklets_plot2(i1,i2,tracklets1(match.subs1),tracklets2(match.subs2));
    
    if i>1
        [cons,pcons] = match.get_cons(pmatch);
        %tracklets_plot2(i1,i2,tracklets1(match.subs1(cons)),tracklets2(match.subs2(cons)));
    end

    [pts1,valid1] = tracklets_tail(tracklets1);
    [pts2,valid2] = tracklets_tail(tracklets2);
    if i==1
        % init the KLT trackers with matched features
        initialize(klt1,pts1,i1);
        initialize(klt2,pts2,i2);
        a{1} = zeros(6,1);
    else
        setPoints(klt1,pts1,valid1);
        setPoints(klt2,pts2,valid2);

        % prepare inputs for error minimization
        X = pmatch.X(:,pcons);
        x = [pts1(match.subs1(cons),:) pts2(match.subs2(cons),:)]';
        v = ~isnan(X(1,:));
        X = X(:,v);
        x = x(:,v);

        % plot inputs
        x0 = h2e(P0*e2h(X));
        figure('units','normalized','outerposition',[0 0 1 1])
        showMatchedFeatures(i1,pi1,x0',x(1:2,:)');
        title('putative matches'); legend('prev','current');
        saveas(gcf,fullfile(DBG_DIR,sprintf('pre_matches_%04d.png',i))); close;
        
        % estimate params
        [a{i},inliers,tr0,predict,rms] = ransac_minimize_reproj(X,x,param);
        
        % plot outputs
        gt = mat2tr(poses(:,:,i));
        x0 = h2e(P0*e2h(X(:,inliers)));
        figure('units','normalized','outerposition',[0 0 1 1])
        showMatchedFeatures(i1,pi1,x0',x(1:2,inliers)','blend');
        title('estimation support set'); legend('prev','current');
        startx = 5; starty = 10; inc = 20;
        text(startx,starty,sprintf('inliers: %d of %d',numel(inliers),size(X,2)),'color','green');
        starty = starty + inc; text(startx,starty,sprintf('residual rms: %g',rms),'Color','green');
        starty = starty + inc; text(startx,starty,sprintf('ground truth param: %s',sprintf('%0.2g ',gt)),'Color','green');
        starty = starty + inc; text(startx,starty,sprintf('current estimation: %s',sprintf('%0.2g ',tinv(a{i}))),'Color','green');
        starty = starty + inc; text(startx,starty,sprintf('initial guess: %s',sprintf('%0.2g ',tinv(tr0))),'Color','green');
        saveas(gcf,fullfile(DBG_DIR,sprintf('post_matches_%04d.png',i))); close;
        
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(211); showMatchedFeatures(i1,i1,predict(1:2,:)',x(1:2,inliers)','blend');
        title('residuals for inliers'); legend('predicted','measured');

        T = tr2mat(a{i});
        % predicted feature locations both inliers & outliers
        pred = h2e(param.K*T(1:3,:)*e2h(X));
        % image measurements in the previous frame
        x0 = h2e(P0*e2h(X));
        im = imfuse(i1,pi1,'blend');
        subplot(212);
        imshow(im); hold on;
        plot(pred(1,:),pred(2,:),'or',x0(1,:),x0(2,:),'+b',x(1,:),x(2,:),'ob');
        % connect predicted and measured
        lx = [pred(1,:); x(1,:)]; ly = [pred(2,:); x(2,:)]; line(lx,ly,'color','yellow');
        % connect putative matches
        lx = [x0(1,:); x(1,:)]; ly = [x0(2,:); x(2,:)]; line(lx,ly,'color','magenta');
        title('residuals for all putative matches'); legend('predicted','prev measured','cur measured');
        saveas(gcf,fullfile(DBG_DIR,sprintf('residuals_%04d.png',i))); close;
    end
    
    %[x,~] = track2cell(tracks);
    %if i>= size(a,1)
    %    sigma = diag(ones(1,4*numel(x)));
    %    [aba,bba] = simple_ba(x,sigma,a,b,param);
    %end
    pi1 = i1; pi2 = i2; pmatch = match.copy();
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

function [valid,patch] = patch_get(im,center,patchr)
% extract a square patch centered @ center of size
% (2*patch_sz+1)x(2*patch_sz+1)
% if the center is too close to the border, valid will be false and
% patch=[]
x1 = center(1)-patchr;
x2 = center(1)+patchr;
y1 = center(2)-patchr;
y2 = center(2)+patchr;
if x1<1 || y1<1 || x2>size(im,2) || y2>size(im,1)
    valid = false;
    patch = [];
else
    patch = double(im(y1:y2,x1:x2));
    valid = true;
end
end

function [patches, coords] = patches_get(im,c,searchr,patchr)

[h,w] = size(im);
is_valid = @(pt) pt(1)-patchr>1 & pt(2)-patchr>1 & pt(1)+patchr<w & pt(2)+patchr<h;
k=1;
for j1=-searchr:searchr
    for j2=-searchr:searchr
        pt = [c(1)+j1;c(2)+j2];
        if is_valid(pt)
            patches(:,k) = double(reshape(im(pt(2)-patchr:pt(2)+patchr,pt(1)-patchr:pt(1)+patchr),[(2*patchr+1)^2 1]));
            coords(:,k) = pt;
            k = k+1;
        end
    end
end

end

function pts = features_pts(features)
valid = find(arrayfun(@(x) x.valid(),features));
len = length(valid);
pts = nan(2,len);
for j=1:len
    pts(:,j) = features(valid(j)).pt;
end

end

function frames = features_init(i1,param)
%corners = corner(i1,param.feat_method,param.feat_num)';
%valids = (corners(1,:)-param.patchr>1) &...
%    (corners(1,:)+param.patchr<size(i1,2)) &...
%    (corners(2,:)-param.patchr>1) &...
%    (corners(2,:)+param.patchr<size(i1,1));
%corners = corners(:,valids);
corners = corner_blob(i1,param.feat_num,param.patchr);

len = size(corners,2);
frames(len) = feat;
for j=1:len
    c = round(corners(:,j));
    [~, patch] = patch_get(i1,c,param.patchr);
    frames(j) = feat(c,patch(:));
end
end

function frames = features_augment(im,pts,param)
% pts is a set of existing features
%corners = corner(im,param.feat_method,param.feat_num)';
%valids = (corners(1,:)-param.patchr>1) &...
%    (corners(1,:)+param.patchr<size(im,2)) &...
%    (corners(2,:)-param.patchr>1) &...
%    (corners(2,:)+param.patchr<size(im,1));
%corners = corners(:,valids);
corners = corner_blob(im,param.feat_num,param.patchr);
valids = false(size(corners,2),1);
for j=1:size(corners,2)
    d = bsxfun(@minus,pts,corners(:,j));
    if all(sum(d.*d)>25)
        valids(j) = true;
    end
end
corners = corners(:,valids);
len = size(corners,2);
frames(len) = feat;
for j=1:len
    c = round(corners(:,j));
    [~, patch] = patch_get(im,c,param.patchr);
    frames(j) = feat(c,patch(:));
end
end

function frames = features_track(im,klt,param)
[corners, validity] = step(klt,im);
len = size(corners,1);
frames(len) = feat;
for j=1:len
    c = round(corners(j,:));
    [valid,patch] = patch_get(im,c,param.patchr);
    if valid==true && validity(j)==true,
        frames(j) = feat(c,patch(:));
    end
end
end


function [X,d,J] = triangulate1(x1,x2,param)
[X, d, J] = triangulate_naive(x1,x2,param.base,param.calib.f,...
    param.calib.cu,param.calib.cv);
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

function tracklets = tracklets_init(features,t)
len = length(features);
tracklets(len) = tracklet;
for j=1:len
    tracklets(j).push_back(features(j),t);
end
end

function tracklets = tracklets_update(tracklets,features,new_features,t)
arrayfun(@(x,y) x.update(y,t),tracklets,features)
tlen = max(size(tracklets));
flen = length(new_features);
tracklets(tlen+flen) = tracklet;
for j=1:flen
    tracklets(tlen+j).push_back(new_features(j),t);
end
end

function [pts,val] = tracklets_tail(tracklets)
pts = [arrayfun(@(x) x.last().x,tracklets); arrayfun(@(x) x.last().y, tracklets)]';
val = arrayfun(@(x) x.valid(),tracklets)';
end

function pts = tracklets_ptail(tracklets)
pts = [arrayfun(@(x) x.plast().x,tracklets); arrayfun(@(x) x.plast().y, tracklets)]';
end

function tracklets_plot(tracklets,i1,i1p,tit)
[pts,valid] = tracklets_tail(tracklets);
ppts = tracklets_ptail(tracklets);
figure; showMatchedFeatures(i1,i1p,pts(valid,:),ppts(valid,:),'blend'); title(tit);
end

function tracklets_plot2(i1,i2,tracklets1,tracklets2,tit)
%% plot matching left/right tracklets.  tracklets1(j) <-> tracklets2(j)
if nargin==4
    tit = '';
end
[pts1,~] = tracklets_tail(tracklets1);
[pts2,~] = tracklets_tail(tracklets2);
figure; showMatchedFeatures(i1,i2,pts1,pts2,'blend'); title(tit);
end

function c = corner_blob(im,feat_num,patchr)
% Geiger StereoScan masks
blob_mask = [-1 -1 -1 -1 -1; -1 +1 +1 +1 -1; -1 +1 +8 +1 -1; -1 +1 +1 +1 -1; -1 -1 0 -1 -1];
corn_mask = [-1 -1 +0 +1 +1; -1 -1  0 +1 +1;  0  0  0  0  0; +1 +1  0 -1 -1; +1 +1 0 -1 -1];

% apply the filters
blobs = imfilter(im,blob_mask);
corns = imfilter(im,corn_mask);

bw = imregionalmax(blobs,8); bw = bwmorph(bw,'shrink',Inf); [blobr,blobc] = find(bw);
bw = imregionalmax(corns,8); bw = bwmorph(bw,'shrink',Inf); [cornr,cornc] = find(bw);

% discard features that are too close to edges
blob_valid = blobr>patchr & blobr<size(im,1)-patchr & blobc>patchr & blobc<size(im,2)-patchr;
corn_valid = cornr>patchr & cornr<size(im,1)-patchr & cornc>patchr & cornc<size(im,2)-patchr;
blobr = blobr(blob_valid); blobc = blobc(blob_valid);
cornr = cornr(corn_valid); cornc = cornc(corn_valid);

% subsample, since we geet too much
if length(blobr)>feat_num
    blob_active = datasample(1:length(blobr),feat_num,'Replace',false);
end
if length(cornr)>feat_num
    corn_active = datasample(1:length(cornr),feat_num,'Replace',false);
end

% format the result
c(1,:) = [blobc(blob_active); cornc(corn_active)]';
c(2,:) = [blobr(corn_active); cornr(corn_active)]';

end

function [features,valid] = tracklets_track(tracklets,i1,searchr,patchr)
% track features using simple template matching
[pts,valid] = tracklets_tail(tracklets);
pts = pts'; valid = valid';
features(length(valid)) = feat;
parfor j=find(valid)
    c = pts(:,j); d = double(tracklets(j).last().descriptor);
    [patches,coords] = patches_get(i1,c,searchr,patchr);
    [~,ind] = min(sum(abs(bsxfun(@minus,patches,d))));
    best_patch = patches(:,ind);
    best_pt = coords(:,ind);
    features(j) = feat(best_pt, best_patch);
end
end