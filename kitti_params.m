function param = kitti_params(P0, P1)

param.P1 = P0;
param.P2 = P1;
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
param.ransac_iter = 500;
% observation dimension
param.obd = 4;
% camera parameter vector dimension
param.ad = 6;
% structure parameter vector dimension
param.bd = 3;
param.threshx = 100;
param.threshy = 2;
param.lm_max_iter = 100;