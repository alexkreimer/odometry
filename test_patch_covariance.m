function test_patch_covariance()

close all;
dbstop if error;

KITTI_HOME = '/home/kreimer/KITTI/dataset';
DBG_DIR = 'debug';

image_dir = fullfile(KITTI_HOME, 'sequences', '00');
poses_file = fullfile(KITTI_HOME, 'poses','00.txt');

first = 1;
last = 100;

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

param.feat_num = 500;                        % number of corners to extract per image
param.feat_method = 'MinimumEigenvalue';     % corner extraction method (not used)
param.patchr = 3;                            % descriptors are of size (2*patchr+1)**2
param.searchr = 100;                         % search window radius for feature tracking

param.min_d = 4;                             % minimal acceptable disparity
param.min_disp = -inf;
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

param.threshx = 100;
param.threshy = 2;

param.lm_max_iter = 1000;
a = repmat({zeros(6,1)},3,1);

% main loop, iteration per image pair
for i=first:last
    % read images
    [i1,i2] = read_kitti_images(image_dir, i);
    
    if i==first
        detector = vis.HarrisCorners(i1,param.feat_num,param.patchr);
        extractor = vis.patch_extractor(i1,detector.corners,param.patchr);
        features = vis.feature(extractor.center,extractor.descriptor);
        tracklets = vis.tracklet(features);
        %tracklets1.plot(i1,sprintf('frame %d: left',i));
    else
        % track features from previous frame into the current frame
        tracklets.track(i1,param.searchr,pi1);
     
        % plot left/right tracks
        %figure; tracklets1.pplot(i1,pi1,sprintf('frame %d: left vs prev',i));

        % add features from the current frame
        tracklets.augment(i1,param.feat_num,param.patchr,5);

        % plot pair matches
        %figure;
        %tracklets_plot2(i1,i2,tracklets1(match.subs1),tracklets2(match.subs2),'stereo pair',false);        
        
        [X,x] = tracklets.get_matched(tracklets2);
        
        % plot inputs
%         x0 = h2e(P0*e2h(X));
%         figure('units','normalized','outerposition',[0 0 1 1])
%         showMatchedFeatures(i1,pi1,x0',x(1:2,:)');
%         title('putative matches'); legend('prev','current');
%         saveas(gcf,fullfile(DBG_DIR,sprintf('pre_matches_%04d.png',i))); close;
        
        % estimate params
        [a{i},inliers,tr0,predict,rms] = ransac_minimize_reproj(X,x,param);
        
        % plot outputs
%         gt = mat2tr(poses(:,:,i));
%         x0 = h2e(P0*e2h(X(:,inliers)));
%         figure('units','normalized','outerposition',[0 0 1 1])
%         showMatchedFeatures(i1,pi1,x0',x(1:2,inliers)','blend');
%         title('estimation support set'); legend('prev','current');
%         startx = 5; starty = 10; inc = 20;
%         text(startx,starty,sprintf('inliers: %d of %d',numel(inliers),size(X,2)),'color','green');
%         starty = starty + inc; text(startx,starty,sprintf('residual rms: %g',rms),'Color','green');
%         starty = starty + inc; text(startx,starty,sprintf('ground truth param: %s',sprintf('%0.2g ',gt)),'Color','green');
%         starty = starty + inc; text(startx,starty,sprintf('current estimation: %s',sprintf('%0.2g ',tinv(a{i}))),'Color','green');
%         starty = starty + inc; text(startx,starty,sprintf('initial guess: %s',sprintf('%0.2g ',tinv(tr0))),'Color','green');
%         saveas(gcf,fullfile(DBG_DIR,sprintf('post_matches_%04d.png',i))); close;
        
%         figure('units','normalized','outerposition',[0 0 1 1])
%         subplot(211); showMatchedFeatures(i1,i1,predict(1:2,:)',x(1:2,inliers)','blend');
%         title('residuals for inliers'); legend('predicted','measured');

        if i == 1
            Tgt = inv(poses_gt(:,:,1));
        else
            T1 = [poses_gt(:,:,i);0 0 0 1];     % represents frame #i as seen from frame #0
            T2 = [poses_gt(:,:,i-1);0 0 0 1];   % represents frame #i-1 as seen from frame #0
            Tgt = inv(T1\T2;
            Tgt = Tgt(1:3,:);
        end
        
        xgt = h2e(param.K*Tgt*e2h(X));
        
        Test = tr2mat(a{i});
        % predicted feature locations both inliers & outliers
        xpred = h2e(param.K*Test(1:3,:)*e2h(X));
        
        % image measurements in the previous frame
        x0 = h2e(P0*e2h(X));
        
        im = imfuse(i1,pi1,'blend');
        imshow(im); hold on;
        plot(xpred(1,:),xpred(2,:),'or',x0(1,:),x0(2,:),'*y',x(1,:),x(2,:),'ob',xgt(1,:),xgt(2,:),'og');
        legend('predicted','initial','measured','gt prediction');
        
        % connect predicted and initial
        lx = [xpred(1,:); x0(1,:)]; ly = [xpred(2,:); x0(2,:)]; line(lx,ly,'color','yellow');

        % connect measured and initial
        lx = [x(1,:); x0(1,:)]; ly = [x(2,:); x0(2,:)]; line(lx,ly,'color','yellow');

        % connect measured and initial
        lx = [xgt(1,:); x0(1,:)]; ly = [xgt(2,:); x0(2,:)]; line(lx,ly,'color','yellow');
        
        title('residuals for all putative matches');
        %saveas(gcf,fullfile(DBG_DIR,sprintf('residuals_%04d.png',i))); close;
        close;
    end
    
    %[x,~] = track2cell(tracks);
    %if i>= size(a,1)
    %    sigma = diag(ones(1,4*numel(x)));
    %    [aba,bba] = simple_ba(x,sigma,a,b,param);
    %end
    pi1 = i1; pi2 = i2;
end

figure; hold on;
plot_poses(poses_gt,zeros(3,1),'b');
estimated = poscell2mat(a);
plot_poses(estimated,zeros(3,1),'r');

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


function tracklets_plot2(i1,i2,tracklets1,tracklets2,titl,separate_figures)
% plot matching left/right tracklets.  tracklets1(j) <-> tracklets2(j)
if nargin<6
    separate_figures = false;
end
if nargin<5
    titl = '';
end
[pts1,~] = tracklets_tail(tracklets1);
[pts2,~] = tracklets_tail(tracklets2);
if separate_figures
    for i=1:size(pts1,2)
        figure;
        showMatchedFeatures(i1,i2,pts1(i,:),pts2(i,:),'montage');
        waitforbuttonpress;
        close;
    end
else
    figure; showMatchedFeatures(i1,i2,pts1,pts2); title(titl);
end
end

function plot_poses(poses,p,ptspec)
t = nan(3,size(poses,3));
for i=1:size(poses,3)
    pos = poses(:,:,i);
    t(1:3,i) = pos*[p;1];
end
plot3(t(1,:),t(2,:),t(3,:),ptspec);
end

function pos = poscell2mat(a)
pos = nan(3,4,length(a));
for i=1:length(a)
    p = tr2mat(a{i});
    pos(:,:,i) = p(1:3,:);
end
end

