% This function estimates the motion of a stereo rig between 2 pairs of stereo frames

function a_est = estimate_stereo_motion(x,tail,K,num_pts,R0,t0,varargin)
% x(1:2,1:num_pts) projections in the left camera, first pose
% x(1:2,num_pts+1:2*num_pts) projections in the right camera, first pose
% x(3:4,1:num_pts) projections if the left camera, 2nd pose
% x(3:4,num_pts+1:2*num_pts) projections in the right camera, 2nd pose
% K intrinsic calibration matrix
% 3d point(s) visible in the first camera, used to resolve the twisted pair
% ambiguity
% R0,t0 extrinsics of the stereo rig

p = inputParser;

addOptional(p,'i1',[]);
addOptional(p,'pi1',[]);
addOptional(p,'pi2',[]);
addOptional(p,'a0',[]);
addOptional(p,'i',[]);
addOptional(p,'DBG_DIR',[]);

parse(p,varargin{:});


% Estimate the motion of the left camera.  We determine all motion parametes
% except the singed magnitue of the translation

x1 = x(1:2,1:num_pts);           % points from i1
x2 = x(3:4,1:num_pts);           % points from pi1

% Estimate 1st essential: x2'*F*x1
%E1 = K'*fundmatrix(x1,x2)*K; 
%[F1(:,:,1), inliers7] = ransacfitfundmatrix7(x1, x2, .001, 1);
[F1, inliers] = estimateFundamentalMatrix(x1', x2', 'NumTrials', 4000);
%[F1(:,:,3), inliers8] = ransacfitfundmatrix(x1,x2,.001,1);

[~,~,v] = svd(F1);
e1 = h2e(v(:,3));

valid = false(length(tail));
for i=1:length(tail)
    if size(tail{i},2)==3
        valid(i) = true;
    end
end
if size([tail{valid}],2)>1
    pts = [tail{valid}];
end

inliers = find(inliers);
r_tot = 0;
for i=1:numel(inliers)
    p2 = [x2(:,inliers(i));1];
    p1 = [x1(:,inliers(i));1];
    r = p2'*F1*p1;
    r = r*r;
    r_tot = r_tot + r;
end

r_rms = sqrt(r/numel(inliers));
fprintf('%d inlier rms (left cam estimate) = %g\n',numel(inliers),r_rms);

E1 = K'*F1*K;

if ~isempty(p.Results.i1) && ~isempty(p.Results.pi1)
    i1 = p.Results.i1;
    pi1 = p.Results.pi1;
    
    figure;
    
    subplot(211); imshow(i1); hold on;
    title('inliers and epipolar lines in the first image');
    plot(x1(1,inliers),x1(2,inliers),'go');
    
    epiLines = epipolarLine(F1', x2(:,inliers)');
    points = lineToBorderPoints(epiLines, size(i1));
    line(points(:, [1,3])', points(:, [2,4])');
    
    subplot(212); imshow(pi1); hold on;
    title('inliers and epipolar lines in the 2nd image');
    plot(x2(1,inliers),x2(2,inliers),'go');
    epiLines = epipolarLine(F1, x1(:,inliers)');
    points = lineToBorderPoints(epiLines, size(pi1));
    line(points(:, [1,3])', points(:, [2,4])');
    
    truesize;
    saveas(gcf,fullfile(p.Results.DBG_DIR,sprintf('epip1_%04d.png',p.Results.i))); close;
end

% Decompose essential
s1 = decompose_essential(E1,K,x(1:4,1:num_pts));

% origin and orientation of the left_cam(1) as seen from left_cam(2)
t1 = s1(1:3,4);
R1 = s1(1:3,1:3);

% if we have ground truth, verify estimation
if ~isempty(p.Results.a0)
    T1 = inv([s1; 0 0 0 1]);
    T  = tr2mat(a0);
    assertT(T1,T);
end

% Estimate motion between initial position of the right camera and current
% position of the left camera
x1 = x(1:2,(num_pts+1):(2*num_pts)); %points from pi2
x2 = x(3:4,1:num_pts); % points from i1

%E2 = K'*fundmatrix(e2h(x1),e2h(x2))*K; % x2'*F*x1
%[F2, inliers] = ransacfitfundmatrix7(x1, x2, .01);
[F2,inliers] = estimateFundamentalMatrix(x1', x2', 'NumTrials', 4000);
E2 = K'*F2*K;

inliers = find(inliers);
r_tot = 0;
for i=1:numel(inliers)
    p2 = [x2(:,inliers(i));1];
    p1 = [x1(:,inliers(i));1];
    
    r = p2'*F2*p1;
    r = r*r;
    r_tot = r_tot + r;
end

r_rms = sqrt(r/numel(inliers));
fprintf('%d inlier rms (right cam estimate) = %g\n',numel(inliers),r_rms);

if ~isempty(p.Results.pi2) && ~isempty(p.Results.i1)
    i1 = p.Results.i1;
    pi2 = p.Results.pi2;
    figure;
    
    subplot(211); imshow(pi2); hold on;
    title('inliers and epipolar lines in the first image');
    plot(x1(1,inliers),x1(2,inliers),'go');
    
    epiLines = epipolarLine(F1', x2(:,inliers)');
    points = lineToBorderPoints(epiLines, size(i1));
    line(points(:, [1,3])', points(:, [2,4])');
    
    subplot(212); imshow(i1); hold on;
    title('inliers and epipolar lines in the 2nd image');
    plot(x2(1,inliers),x2(2,inliers),'go');
    epiLines = epipolarLine(F1, x1(:,inliers)');
    points = lineToBorderPoints(epiLines, size(pi2));
    line(points(:, [1,3])', points(:, [2,4])');
    truesize;
    saveas(gcf,fullfile(p.Results.DBG_DIR,sprintf('epip2_%04d.png',p.Results.i))); close;
end

s2 = decompose_essential(E2,K,[x1;x2]);

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

t2 = s2(1:3,4);
R2 = s2(1:3,1:3);

t2 = -R2'*t2;     % t2 in the frame of the left camera
t1 = -R1'*t1;

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
T1 = [R1' t1;0 0 0 1];

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
