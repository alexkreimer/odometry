% This function estimates the motion of a stereo rig between 2 pairs of stereo frames

function a_est = estimate_stereo_motion(x,K,num_pts,R0,t0,varargin)
% x(1:2,1:num_pts) feature coords as seen in the left camera, pose {i-1}
% x(3:4,1:num_pts) same features coords the left camera, pose {i}
%
% x(1:2,num_pts+1:2*num_pts) feature coords in the right camera, {i-1}
% x(3:4,num_pts+1:2*num_pts) feature coords in the right camera, {i}
%
% K camera intrinsics
% R0,t0 stereo rig extrinsics
%
% a_est estimated 6-dof motion parameter vector. Describes pose {i-1} as seen
% in frame {i}

p = inputParser;
addOptional(p,'i1',[]);
addOptional(p,'pi1',[]);
addOptional(p,'pi2',[]);
addOptional(p,'a0',[]);
addOptional(p,'i',[]);
addOptional(p,'DBG_DIR',[]);
addOptional(p,'do_cr', false, @islogical);
addOptional(p,'gt',[])
addOptional(p,'ind',[])

parse(p,varargin{:});

% Estimate the motion of the left camera.  We determine all motion parametes
% except the singed magnitue of the translation

x1 = x(3:4, 1:num_pts); % points from prev i1
x2 = x(1:2, 1:num_pts); % points from i1

[F1_est, ~, T1_est, inliers1] = estimateF(x1, x2, K);  % estimate E/F and decompose the essential
if ~isempty(p.Results.gt) && ~isempty(p.Results.i1) && ~isempty(p.Results.pi1)
    dbg_estimation(F1_est, K, T1_est, inv([p.Results.gt; 0 0 0 1]), x1, x2, inliers1, p.Results.pi1, p.Results.i1,...
        'prev left', 'left', p.Results.DBG_DIR, p.Results.ind, 1, 0);
end

% if we have ground truth, verify estimation
if ~isempty(p.Results.a0)
    T  = tr2mat(a0);
    assertT(T_est,T);
end

% Estimate motion between initial position of the right camera and current
% position of the left camera
x1 = x(3:4,(num_pts+1):(2*num_pts)); % points from pi2
x2 = x(1:2,1:num_pts);               % points from i1
[F2_est, ~, T2_est, inliers2] = estimateF(x1, x2, K);  % estimate E/F and decompose the essential
if ~isempty(p.Results.gt) && ~isempty(p.Results.i1) && ~isempty(p.Results.pi2)
    T2_gt = [p.Results.gt; 0 0 0 1];
    T0 = [R0,t0; 0 0 0 1];
    T2_gt = inv(T2_gt*T0);
    dbg_estimation(F2_est, K, T2_est, T2_gt, x1, x2, inliers2, p.Results.pi2, p.Results.i1, 'prev right', 'current left', p.Results.DBG_DIR, p.Results.ind, 2, 1);
end

if ~isempty(p.Results.a0)
    % pose of lcam2 in rcam1
    T2 = [s2; 0 0 0 1];
    
    % rcam1 described in lcam1
    T0 = [R0,t0;0 0 0 1];
    
    % lcam2 in rcam1
    T = T0\T;
    T = unit_t(T);
    
    assertT(T2,T,'inv');
end

t1 = T1_est(1:3,4);
R1 = T1_est(1:3, 1:3);
t2 = T2_gt(1:3,4);
R2 = T2_gt(1:3,1:3);

%t2 = -R2'*t2;     % t2 in the frame of the left camera
%t1 = -R1'*t1;

% now solve for scale of the motion
A = [t1, -t2];
b = t0;
[U,S,V] = svd(A);
b = U'*b;
d = diag(S);
y = nan(size(A,2),1);
for j=1:length(d)
    y(j) = b(j)/d(j);
end
c = V*y;

% form the output
t1 = t1*c(1);
T1 = [R1 t1;0 0 0 1];

figure;
quiver3(0,0,0,t1(1),t1(2),t1(3)); hold on;
quiver3(t0(1),t0(2),t0(3),t2(1),t2(2),t2(3));

% if we have ground truth, verify estimation
if ~isempty(p.Results.a0)
    % ground truth is a pose of cam2 in cam1
    assertT(T1,tr2mat(a0));
end

a_est = mat2tr(T1);

end

function T = unit_t(T)
T(1:3,4) = T(1:3,4)/norm(T(1:3,4));
end

function err = rot_error(pose_error)
a = pose_error(1,1);
b = pose_error(2,2);
c = pose_error(3,3);
d = 0.5*(a+b+c-1.0);
err = acos(max(min(d,1.0),-1.0));
end

function err = trans_error(pose_error)
dx = pose_error(1,4);
dy = pose_error(2,4);
dz = pose_error(3,4);

err = sqrt(dx*dx+dy*dy+dz*dz);
end

function plot_epip(F, x1, x2, i1, i2, title1, title2)
figure;
subplot(211); imshow(i1); hold on;
title(title1);
plot(x1(1,:), x1(2,:),'go');
epiLines = epipolarLine(F', x2');
points = lineToBorderPoints(epiLines, size(i1));
line(points(:, [1,3])', points(:, [2,4])');
subplot(212); imshow(i2); hold on;
title(title2);
plot(x2(1,:), x2(2,:), 'go');
epiLines = epipolarLine(F, x1');
points = lineToBorderPoints(epiLines, size(i2));
line(points(:, [1,3])', points(:, [2,4])');
truesize;
end

function err = residual_error(F, x1, x2, inliers)
% computes average error from putative matches to the correspnding epipolar
% lines

if nargin==3
    inliers = 1:length(x1);
end

err=0;
for i = inliers(:)'
    l2 = F*[x1(:, i); 1];
    l2 = l2/sqrt(l2(1)*l2(1)+l2(2)*l2(2));
    d1 = [x2(:, i); 1]'*l2;
    
    l1 = F'*[x2(:, i); 1];
    l1 = l1/sqrt(l1(1)*l1(1)+l1(2)*l1(2));
    
    d2 = [x1(:, i); 1]'*l1;
    err = err + d1+d2;
end

err = err/numel(inliers);
end


function dbg_estimation(F_est, K, T_est, T_gt, x1, x2, inliers, i1, i2, title1, title2, DBG_DIR, i, j, single_save)

t_gt = T_gt(1:3, 4);
R_gt = T_gt(1:3, 1:3);
F_gt = inv(K')*skew(t_gt)*R_gt*inv(K);
err = residual_error(F_est, x1, x2);
fprintf('GT: left vs prev left: average residual error %g [px] for %d matches\n', err, length(x1));

pose_error = T_est\T_gt;
fprintf('orientation error: %g [rad]\n', rot_error(pose_error));

plot_epip(F_est, x1(:, inliers), x2(:, inliers), i1, i2, title1, title2);
saveas(gcf,fullfile(DBG_DIR, sprintf('epip_%04d_est_%d.png', i, j))); close;
plot_epip(F_gt, x1(:, inliers), x2(:, inliers), i1, i2, sprintf('GT: %s', title1),...
    sprintf('GT: %s', title2));
saveas(gcf,fullfile(DBG_DIR, sprintf('epip_%04d_gt_%d.png', i, j))); close

if single_save
    for k = inliers(:)'
        figure; ax = axes;        
        showMatchedFeatures(i1, i2, x1(:, k)', x2(:, k)', 'montage', 'Parent', ax);
        title(ax, 'Candidate point matches');
        saveas(gcf, fullfile(DBG_DIR, sprintf('matches_%04d_%d_%d.png', i, j, k)));
        close
    end
else
    figure; ax = axes;
    showMatchedFeatures(i1, i2, x1(:,inliers)', x2(:,inliers)', 'Parent', ax);
    title(ax, 'Candidate point matches'); legend(ax, 'Matched points 1','Matched points 2');
    saveas(gcf,fullfile(DBG_DIR, sprintf('matches_%04d_%d.png', i, j))); close
end

end

function [F_est, E_est, T_est, inliers] = estimateF(x1, x2, K, thresh)

if nargin==3
    thresh = 2;
end
% Estimate fundamental: x2'*F*x1
tic; [F_est, inliers] = estimateFundamentalMatrix(x1', x2', 'DistanceType', 'Sampson', 'Method', 'RANSAC', 'DistanceThreshold', thresh, 'NumTrials', 4000); toc;
inliers = find(inliers);
err = residual_error(F_est, x1, x2, inliers);
fprintf('fundamental estimation, left vs prev left: average residual error %g [px] for %d inliers\n', err, numel(inliers));
% compute essential
E_est = K'*F_est*K;

% Decompose essential
T_est = decompose_essential(E_est, K, [x1;x2]);
T_est = inv([T_est;0 0 0 1]);
end
