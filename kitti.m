function kitti()

close all;
dbstop if error;

KITTI_HOME = '/home/kreimer/KITTI/dataset';
DBG_DIR = 'debug';

image_dir = fullfile(KITTI_HOME, 'sequences', '00');
poses_file = fullfile(KITTI_HOME, 'poses','00.txt');

first = 1;
last  = 3;

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
a = repmat({zeros(6,1)},last,1);
a_est = repmat({zeros(6,1)},last,1);
% main loop, iteration per image pair
for i=first:last
    % read images
    [i1,i2] = read_kitti_images(image_dir, i);
    pos_3d(:,:,1) = tr2mat(zeros(6,1));
    poses(:,:,1) = tr2mat(zeros(6,1));

    if i==first
        detector1 = vis.HarrisCorners(i1,param.feat_num,param.patchr);
        detector2 = vis.HarrisCorners(i2,param.feat_num,param.patchr);
        
        extractor1 = vis.patch_extractor(i1,detector1.corners,param.patchr);
        extractor2 = vis.patch_extractor(i2,detector2.corners,param.patchr);
        
        features1 = vis.feature(extractor1.center,extractor1.descriptor);
        features2 = vis.feature(extractor2.center,extractor2.descriptor);
        
        tracklets1 = vis.tracklet(features1);
        tracklets2 = vis.tracklet(features2);
        
%        tracklets1.plot(i1,sprintf('frame %d: left',i));
%        tracklets2.plot(i2,sprintf('frame %d: right',i));

        % match between current left/right pair
        tracklets1.hmatch(tracklets2,param);
else
        % track features from previous frame into the current frame
        % these are intra-view tracks (no calibration available)
        tracklets1.track(i1,param.searchr,pi1);
        tracklets2.track(i2,param.searchr,pi2);
        
        % plot left/right tracks
        %        figure; 
        %        subplot(211); tracklets1.pplot(i1,pi1,sprintf('frame %d: left vs prev',i));
        %        subplot(212); tracklets2.pplot(i2,pi2,sprintf('frame %d: right vs prev',i));
%        saveas(gcf,fullfile(DBG_DIR,sprintf('feature_tracks_%04d.png',i))); close;
        
        % add features from the current frame
        tracklets1.augment(i1,param.feat_num,param.patchr,5);
        tracklets2.augment(i2,param.feat_num,param.patchr,5);
        
        % match between current left/right pair
        tracklets1.hmatch(tracklets2,param);
        
%        figure; tracklets1.plot_circles(i1,pi1,i2,pi2,tracklets2);
%        saveas(gcf,fullfile(DBG_DIR,sprintf('feature_circles_%04d.png',i))); close;

%        if i>2
%            tracklets1.fit_line();
%        end
        
        % plot pair matches
        %figure;
        %tracklets_plot2(i1,i2,tracklets1(match.subs1),tracklets2(match.subs2),'stereo pair',false);
        
        %[Xp,xp,x,ind,tail] = tracklets1.get_matched(tracklets2);
        
        % estimate params
        %[a{i},inliers,tr0,predict,rms] = ransac_minimize_reproj(Xp,x,param);
        %if i>1
        %    pos_3d(:,:,i) = inv(tr2mat(a{i})*pos_3d(:,:,i-1));
        %end
        %eval_estimate(poses_gt(:,:,i),pos_3d(:,:,i),Xp,xp,x,inliers,tr0,predict,rms,pi1,i1,pi2,i2,param.K,i,DBG_DIR);
        
        %num_pts = size(x,2);
        %x = [x(1:2,:) x(3:4,:); xp(1:2,:) xp(3:4,:)];
        %a_est{i} = estimate_stereo_motion(x,param.K,num_pts,eye(3),[param.base,0,0]','i1',i1,'pi1',pi1,'pi2',pi2,'i',i,'DBG_DIR',DBG_DIR);
        %a_est{i} = estimate_stereo_motion(x,param.K,num_pts,eye(3),[param.base,0,0]','DBG_DIR',DBG_DIR);
        %if i>1
        %    poses(:,:,i) = inv(tr2mat(a_est{i})*poses(:,:,i-1));
        %end
    end
 
    pi1 = i1; pi2 = i2;
end

save('tracklets','tracklets1', 'tracklets2');
savePoses('poses.txt', poses);
%save('patches.mat','data');
estimated = poscell2mat(a);
estimated1 = poscell2mat(a_est);

N = size(estimated,3);

figure; hold on;
plot_poses(poses_gt,N,'b');


plot_poses(estimated,N,'r');


plot_poses(estimated1,N,'g');

legend('gt','3d est','geom est');
end

function [i1,i2] = read_kitti_images(seq_home, idx)
file_names = {fullfile(seq_home, 'image_0', sprintf('%06d.png',idx))
    fullfile(seq_home, 'image_1', sprintf('%06d.png',idx))};
i1 = imread(file_names{1});
i2 = imread(file_names{2});
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

function plot_poses(poses,n,ptspec)
t = nan(3,size(poses,3));
for i=1:min(n,size(poses,3))
    if size(poses,1) == 3
        pos = [poses(:,:,i);0,0,0,1];
    else
        pos = poses(:,:,i);
    end
    t(1:3,i) = h2e(pos*[0;0;0;1]);
end
plot(t(1,:),t(3,:),ptspec);
plot(t(1,:),t(3,:),'og');
end

function pos = poscell2mat(a)
pos = nan(3,4,length(a));
for i=1:length(a)
    p = tr2mat(a{i});
    if i>1
        p = p*[pos(:,:,i-1);0,0,0,1];
        pos(:,:,i) = p(1:3,:);
    else
        pos(:,:,i) = p(1:3,:);
    end
end

for i=1:length(a)
    p = inv([pos(:,:,i);0 0 0 1]);
    pos(:,:,i) = p(1:3,:);
end

end

function [r_err,t_err] = eval_estimate(pos_gt,pos_est,Xp,xp,x,inliers,tr0,predict,rms,pi1,i1,pi2,i2,K,i,DBG_DIR)

pos_gt  = [pos_gt;0 0 0 1];
pos_err = mat2tr(inv(pos_est)*pos_gt);

r_err = acos(min(1.0,sum(abs(pos_err(1:3)))));
t_err = norm(pos_err(4:6))/norm(pos_gt(1:3,4));

figure; 
subplot(211);
hold on;
title(sprintf('measurement vs. prediction, rms=%g',rms));
imshow(i1);
hold on;
plot(x(1,inliers),x(2,inliers),'or');
plot(predict(1,:),predict(2,:),'.g');
legend('measurement','prediction');
subplot(212);
imshow(i2);
hold on;
plot(x(3,inliers),x(4,inliers),'or');
plot(predict(3,:),predict(4,:),'.g');
legend('measurement','prediction');
saveas(gcf,fullfile(DBG_DIR,sprintf('estimation_%04d.png',i))); close;

end
