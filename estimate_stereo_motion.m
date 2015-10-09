% This function estimates the motion of a stereo rig between 2 pairs of stereo frames

function a_est = estimate_stereo_motion(x, K, num_pts, R0, t0, param, varargin)
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

debug = 1;

fprintf('1st estimation\n');
if debug
    T_gt = [p.Results.gt; 0 0 0 1];

    X = nan(4, num_pts);
    for i = 1:num_pts
         X(:,i) = vgg_X_from_xP_nonlin([x(3:4, i) x(3:4, num_pts+i)], {param.P1, param.P2});
    end
    X = h2e(X);
    observe = [x(1:2, 1:num_pts); x(1:2, (num_pts+1):end)];
    [a_est, ~, ~, ~, ~] = ransac_minimize_reproj(X, observe, param);
    T_est1 = tr2mat(a_est);
    E_est1 = skew(T_est1(1:3,4))*T_est1(1:3,1:3);
    F_est1 = inv(K')*E_est1*inv(K);
    
    fig = figure;
    imshow(repmat(p.Results.pi1, [1 1 3]));
    hold on;
    scatter(x(3,1:num_pts), x(4, 1:num_pts), [], X(3,:));
    title('features colored according to their depth, all values le 0 map to 0, ge 50 map to 50');
    caxis([0 50]);
    colormap jet;
    hold off;
    dcm_obj = datacursormode(fig);
    myupdatefcn1 = @(emp, event_obj) myupdatefcn(emp, event_obj, x(:,1:num_pts), X);
    set(dcm_obj, 'UpdateFcn', myupdatefcn1)

    param1 = struct('i1', p.Results.pi1, 'i2', p.Results.i1, 'i1_name',...
        'prev left', 'i2_name', 'left', 'T_gt', T_gt, 'dbg_save', 0,...
        'single_save', 0, 'DBG_DIR', p.Results.DBG_DIR, 'ind1', p.Results.ind,...
        'ind2', 1, 'F_est1', F_est1, 'E_est1', E_est1);
    
    [F1_est, ~, T1_est, inliers1] = estimateF(x1, x2, K, 2, param1);  % estimate E/F and decompose the essential
else
    [F1_est, ~, T1_est, inliers1] = estimateF(x1, x2, K, 2);  % estimate E/F and decompose the essential
end

% if we have ground truth, verify estimation
if ~isempty(p.Results.a0)
    T  = tr2mat(a0);
    assertT(T_est,T);
end

fprintf('2st estimation\n');

% Estimate motion between initial position of the right camera and current
% position of the left camera
x1 = x(3:4,(num_pts+1):(2*num_pts)); % points from pi2
x2 = x(1:2,1:num_pts);               % points from i1

if debug
    T2_gt = [p.Results.gt; 0 0 0 1];
    T0 = [R0,t0; 0 0 0 1];
    T2_gt = inv(T2_gt*T0);
    
    param2 = struct('i1', p.Results.pi2, 'i2', p.Results.i1, 'i1_name',...
        'prev right', 'i2_name', 'left', 'T_gt', T2_gt, 'dbg_save', 0,...
        'single_save', 1, 'DBG_DIR', p.Results.DBG_DIR, 'ind1', p.Results.ind,...
        'ind2', 1, 'F_est1', F_est1, 'E_est1', E_est1);
    [F2_est, ~, T2_est, inliers2] = estimateF(x1, x2, K, 2, param2);  % estimate E/F and decompose the essential
else
    [F2_est, ~, T2_est, inliers2] = estimateF(x1, x2, K, 2);  % estimate E/F and decompose the essential
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
t2 = T2_est(1:3,4);
R2 = T2_est(1:3,1:3);

% t2 = -R2'*t2;     % t2 in the frame of the left camera
% t1 = -R1'*t1;

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
epiLines = epipolarLine(F, x2');
points = lineToBorderPoints(epiLines, size(i1));
line(points(:, [1,3])', points(:, [2,4])');
subplot(212); imshow(i2); hold on;
title(title2);
plot(x2(1,:), x2(2,:), 'go');
epiLines = epipolarLine(F', x1');
points = lineToBorderPoints(epiLines, size(i2));
line(points(:, [1,3])', points(:, [2,4])');
%truesize;
end

function err = residual_error(F, x1, x2, inliers, param, varargin)
% computes average error from putative matches to the correspnding epipolar
% lines

p = inputParser;
addOptional(p, 'dist_fn', @dist_sampson);
parse(p,varargin{:});

if nargin==3
    inliers = 1:length(x1);
end

dist_fn = p.Results.dist_fn;

d = dist_fn(F, x1(:, inliers), x2(:, inliers));
err = sum(abs(d(:)))/length(inliers);

if nargin<5
    return
end

figure;
width1 = .5;
[N, edges] = histcounts(d(1,:));
bar(edges(1:end-1), N, width1, 'FaceColor', [0.2,0.2,0.5]);

if size(d,1) > 1
    hold on;
    [N, edges] = histcounts(d(2,:));
    bar(edges(1:end-1), N, width1/2, 'FaceColor',[0,0.7,0.7],'EdgeColor',[0,0.7,0.7]);
    legend(param.i1_name, param.i2_name);
    title(sprintf('signed distance to the epipolar lines, average %f', err));
else
    legend(param.i1_name);
    title(sprintf('sampson distance, average %f', err));
end

hold off;

end

% symmetric signed distance from pts to the epilines
function d = dist_symm_epiline(F, x1, x2)

n1 = size(x1, 2);
n2 = size(x2, 2);
assert(n1 == n2);
d = nan(2, n1);
for i = 1:n1
    l2 = F*[x1(:, i); 1];
    l2 = l2/sqrt(l2(1)*l2(1)+l2(2)*l2(2));
    d1 = [x2(:, i); 1]'*l2;
    d(1, i) = d1;
    
    l1 = F'*[x2(:, i); 1];
    l1 = l1/sqrt(l1(1)*l1(1)+l1(2)*l1(2));
    
    d2 = [x1(:, i); 1]'*l1;
    d(2, i) = d2;
end

end

function d = dist_sampson(F, x1, x2)
pts1h = e2h(x1);
pts2h = e2h(x2);

pfp = (pts2h' * F)';
pfp = pfp .* pts1h;
d = sum(pfp, 1) .^ 2;

epl1 = F * pts1h;
epl2 = F' * pts2h;
d = d ./ (epl1(1,:).^2 + epl1(2,:).^2 + epl2(1,:).^2 + epl2(2,:).^2);

end

function [F_est, E_est, T_est, inliers] = estimateF(x1, x2, K, thresh, param1)

if nargin == 3
    thresh = 2; %px
end
% Estimate fundamental: x2'*F*x1
[F_est, inliers] = estimateFundamentalMatrix(x1', x2', 'DistanceType',...
    'Sampson', 'Method', 'RANSAC', 'DistanceThreshold', thresh, 'NumTrials', 4000);
inliers = find(inliers);

% compute essential
E_est = K'*F_est*K;

% Decompose essential
T_est = decompose_essential(E_est, K, [x1;x2]);

if nargin<5
    return
end

% use 5-point alg as a check
% x1_norm = nan(3, length(x1));
% x2_norm = nan(3, length(x2));
% for i = 1:size(x1,2)
%     x1_norm(:, i) = K\[x1(:, i); 1];
%     x2_norm(:, i) = K\[x2(:, i); 1];
% end
% evec = calibrated_fivepoint(x1_norm, x2_norm);

F_est1 = param1.F_est1;
E_est1 = param1.E_est1;

% note that the objective optimizes sampson distance and the error here
% below is a epipolar distance, which means that the numeric values will
% probably be different
err = residual_error(F_est, x1, x2, inliers, param1);
fprintf('F_est %s vs %s: residual error %g [px] for %d inliers\n',...
    param1.i1_name, param1.i2_name, err, numel(inliers));

% reconstruct F back from R,t
E_rec = skew(T_est(1:3,4))*T_est(1:3,1:3);
F_rec = inv(K')*E_rec*inv(K);

% calc residuals for the reconstructed F
err = residual_error(F_rec, x1, x2, inliers, param1);
fprintf('F_rec %s vs %s: residual error %g [px] for %d inliers\n',...
    param1.i1_name, param1.i2_name, err, numel(inliers));

% eigenvalue constraint of the essential says that it has to have 2
% identical eigenvalues
lambda = eig(E_est); lambda = lambda/lambda(2);
e_est = null(F_est); e_est = e_est/e_est(3);
fprintf('epipole of F_est: (%f, %f), eigen values of E_est: %f, %f, %f\n',...
    e_est(1), e_est(2), lambda(1), lambda(2), lambda(3));

lambda = eig(E_rec); lambda = lambda/lambda(2);
e_rec = null(F_rec); e_rec = e_rec/e_rec(3);
fprintf(['epipole of F_rec: (%f,%f) eigenvalues of E_rec (essential ',...
    'reconstructed after decomposition): %f, %f, %f\n'], e_rec(1), e_rec(2),...
    lambda(1), lambda(2), lambda(3));

t_gt = param1.T_gt(1:3, 4);
R_gt = param1.T_gt(1:3, 1:3);
E_gt = skew(t_gt)*R_gt;
F_gt = inv(K')*E_gt*inv(K);
err = residual_error(F_gt, x1, x2, inliers, param1);
fprintf('F_gt: %s vs %s: residual error %g [px] for %d matches\n',...
    param1.i1_name, param1.i2_name, err, length(inliers));

lambda = eig(E_gt); lambda = lambda/lambda(2);
e_gt = null(F_gt); e_gt = e_gt/e_gt(3);
fprintf('epipole of F_gt: (%f %f) eigen values of E_gt: %f, %f, %f\n',...
    e_gt(1), e_gt(2), lambda(1), lambda(2), lambda(3));

err = residual_error(F_est1, x1, x2, inliers, param1);
fprintf('F_est1 %s vs %s: residual error %g [px] for %d inliers\n',...
    param1.i1_name, param1.i2_name, err, numel(inliers));

% eigenvalue constraint of the essential says that it has to have 2
% identical eigenvalues
lambda = eig(E_est1); lambda = lambda/lambda(2);
e_est1 = null(F_est1); e_est1 = e_est1/e_est1(3);
fprintf('epipole of F_es1t: (%f, %f), eigen values of E_est: %f, %f, %f\n',...
    e_est1(1), e_est1(2), lambda(1), lambda(2), lambda(3));

fprintf('epipole distance e_rec, e_est: %f\n', norm(e_rec - e_est));
fprintf('epipole distance e_gt, e_est: %f\n', norm(e_gt - e_est));
fprintf('epipole distance e_gt, e_est1: %f\n', norm(e_est1 - e_gt));
fprintf('epipole distance e_est, e_est1: %f\n', norm(e_est1 - e_est));

T_est = inv([T_est;0 0 0 1]);
pose_error = T_est\param1.T_gt;
fprintf('orientation error: %g [rad]\n', rot_error(pose_error));

plot_epip(F_est, x1(:, inliers), x2(:, inliers), param1.i1, param1.i2,...
   sprintf('est: %s', param1.i1_name), sprintf('est: %s', param1.i2_name));

if param1.dbg_save
    save_dbg(fullfile(param1.DBG_DIR, sprintf('epip_%04d_est_%d.png', param1.ind1,...
        param1.ind2)));
    close;
end

plot_epip(F_gt, x1(:, inliers), x2(:, inliers), param1.i1, param1.i2,...
    sprintf('GT: %s', param1.i1_name), sprintf('GT: %s', param1.i2_name));

if param1.dbg_save
    save_dbg(fullfile(param1.DBG_DIR, sprintf('epip_%04d_gt_%d.png', param1.ind1,...
        param1.ind2)));
    close
end

if param1.single_save
    for k = inliers(:)'
        figure; ax = axes;        
        showMatchedFeatures(param1.i1, param1.i2, x1(:, k)', x2(:, k)', 'montage', 'Parent', ax,...
            'PlotOptions', {'ro','go','y-'});
        title(ax, 'Candidate point matches');
        file_name = fullfile(param1.DBG_DIR, sprintf('matches_%04d_%d_%d.png', param1.ind1, param1.ind2, k));
        save_dbg(file_name);
        close;
    end
else
    figure; ax = axes;
    showMatchedFeatures(param1.i1, param1.i2, x1(:,inliers)', x2(:,inliers)', 'Parent',ax);
    hold on;
    plot(e_est(1), e_est(2), 'g*', e_est1(1), e_est1(2), 'r*', e_gt(1), e_gt(2), 'b*');
    legend('matches1', 'matches2', 'estimated', 'estimated1', 'gt');
    title(ax, 'Candidate point matches');
    if param1.dbg_save
        save_dbg(fullfile(DBG_DIR, sprintf('matches_%04d_%d.png', param1.ind1, param1.ind2)));
        close
    end
end

end

function txt = myupdatefcn(empt, event_obj, x, X)
% Customizes text of data tips

pos = get(event_obj,'Position');
d = x(3:4,:) - repmat(pos',[1 length(x)]);
[~, i] = min(sum(d.*d));
txt = {['Depth: ', num2str(X(3,i))], ['x: ', sprintf('%d %d', x(1,i), x(3,i))], ['d: ', sprintf('x:%d y:%d', x(1,i)-x(3,i), x(2,i)-x(4,i))]};
end
