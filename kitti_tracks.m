function kitti_tracks()

close all;
dbstop if error;

KITTI_HOME = '/home/kreimer/KITTI/dataset';
KITTI_HOME = fullfile('F:', 'KITTI' , 'dataset');
DBG_DIR = 'debug';

image_dir  = fullfile(KITTI_HOME, 'sequences', '00');
poses_file = fullfile(KITTI_HOME, 'poses','00.txt');

first = 1;
last  = 20;

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

param.feat_num = 2000;                       % number of corners to extract per image
param.feat_method = 'MinimumEigenvalue';     % corner extraction method (not used)
param.patchr = 3;                            % descriptors are of size (2*patchr+1)**2
param.searchr = 70;                          % search window radius for feature tracking
param.search_thresh = .9;                    % normalized cross correlation threshold, applied in feature tracking

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

param.do_dbg = true;

% main loop, iteration per image pair
for i=first:last
    % read images
    [i1,i2] = read_kitti_images(image_dir, i);
    fprintf('step %d\n', i);
    if i==first
        % first frame, create detectors/descriptors
        detector1 = vis.HarrisCorners(i1,param.feat_num,param.patchr);
        detector2 = vis.HarrisCorners(i2,param.feat_num,param.patchr);
        
        extractor1 = vis.patch_extractor(i1,detector1.corners,param.patchr);
        extractor2 = vis.patch_extractor(i2,detector2.corners,param.patchr);
        
        % extract features
        features1 = vis.feature(extractor1.center,extractor1.descriptor);
        features2 = vis.feature(extractor2.center,extractor2.descriptor);
        
        % init tracks with freshly detected features
        tracklets1 = vis.tracklet(features1);
        tracklets2 = vis.tracklet(features2);

        if param.do_dbg
            figure;
            subplot(211); tracklets1.plot(i1, sprintf('features init, left frame %d',i));
            subplot(212); tracklets2.plot(i2,sprintf('features init, right frame %d',i));
            filename = fullfile(DBG_DIR, sprintf('new_features_%04d.png',i));
            save_dbg(filename);
            close;
        end

        % match between current left/right pair
        tracklets1.hmatch(tracklets2, param);
        
        if param.do_dbg
            tracklets1.plot_matches(tracklets2, i1, i2, sprintf('left vs. right %d', i));
            filename = fullfile(DBG_DIR, sprintf('epipolar_matches_%04d.png',i));
            save_dbg(filename);
            close;
        end
else
        % track features from previous frame into the current frame
        % these are intra-view tracks (no calibration available)
        tracklets1.track(i1, param, pi1);
        tracklets2.track(i2, param, pi2);
        
        if param.do_dbg
            % plot left/right tracks
            tracklets1.pplot1(i1, pi1, sprintf('frame %d: left',i), sprintf('frame %d: prev left',i));
            filename = fullfile(DBG_DIR,sprintf('left_tracks_%04d.png',i));
            save_dbg(filename);
            close;            
            tracklets2.pplot1(i2, pi2, sprintf('frame %d: right',i), sprintf('frame %d: prev right',i));
            filename = fullfile(DBG_DIR,sprintf('right_tracks_%04d.png',i));
            save_dbg(filename);
            close;
        end
        
        % match between current left/right pair
        tracklets1.hmatch(tracklets2,param);

        if param.do_dbg
            tracklets1.plot_matches(tracklets2, i1, i2, sprintf('frame %d: match of left/right under epipolar constraint', i));
            filename = fullfile(DBG_DIR, sprintf('epipolar_matches_%04d.png',i));
            save_dbg(filename);
            close;
        end

        % add features from the current frame
        tracklets1 = tracklets1.augment(i1, param.feat_num,param.patchr, 5);
        tracklets2 = tracklets2.augment(i2, param.feat_num,param.patchr, 5);
        
        if param.do_dbg
            figure; k = tracklets1.plot_circles(i1,pi1,i2,pi2,tracklets2);
            save_dbg(fullfile(DBG_DIR,sprintf('circles_%04d.png',i)));
            close;
            fprintf('%d circles\n', k);
            
            figure; k = tracklets1.plot_circles1(i1,pi1,i2,pi2,tracklets2);
            filename = fullfile(DBG_DIR,sprintf('circles1_%04d.png',i));
            save_dbg(filename);
            close;
        end

        % plot pair matches
        % figure;
        % tracklets_plot2(i1,i2,tracklets1(match.subs1),tracklets2(match.subs2),'stereo pair',false);
    end
    pi1 = i1; pi2 = i2;

end
    save('tracklets','tracklets1', 'tracklets2');
end

