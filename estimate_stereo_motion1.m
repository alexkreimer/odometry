% This function estimates the motion of a stereo rig between 2 pairs of stereo frames

function [a_est, pout] = estimate_stereo_motion1(pin)
% param.c1p feature seen in the left camera, time {i-1}
% param.c1  same features, same camera, time {i}
%
% param.c2p feature coords in the right camera, {i-1}
% param.c2 same features, same camera, time {i}
%
% param.K camera intrinsics
% param.R0, param.t0 stereo rig extrinsics
%
% a_est estimated 6-dof motion parameter vector. Describes pose {i-1} as seen
% in frame {i}

% Left camera motion estimation.
% Determine all motion parametes except the singed magnitue of the translation
[F1, T1, E1, inliers1] = estimateF(pin.c1p, pin.c1, pin.K, 2);  % estimate E/F and decompose the essential

pout.est1.F = F1;
pout.est1.T = T1;
pout.est1.E = E1;
pout.est1.inliers = inliers1;
pout.est1.x1 = pin.c1p;
pout.est1.x2 = pin.c1;
pout.est1.name = 'previous left vs. current left';
pout.xr = 0;

% compute cross-ratios
e = null(F1);
e = h2e(e);
n = size(pin.tracksx,2);
xr = nan(n,1);
yr = nan(n,1);
for i = 1:n
    vx = e(1);
    x  = pin.tracksx(1,i);
    xtt= pin.tracksx(2,i);
    xt = pin.tracksx(3,i);
    xr(i) = (xt-x)*(xtt-vx)/((xtt-xt)*(x-vx));
    
    vy = e(1);
    y  = pin.tracksx(1,i);
    ytt= pin.tracksx(2,i);
    yt = pin.tracksx(3,i);
    yr(i) = (yt-y)*(ytt-vy)/((ytt-yt)*(y-vy));
end

if n>0
    xr(isinf(xr)) = [];
    yr(isinf(yr)) = [];

    % use both x corrdinate and y corrdinate ratios, same same
    data = [xr;yr];
    pd = fitdist(data,'Normal');
    
    % plot the distribution
    x_values = min(data):.1:max(data);
    y = pdf(pd, x_values);
    plot(x_values,y,'LineWidth',2)

    % save the result for output
    pout.xr = pd.mu;
end


% Estimate motion between initial position of the right camera and current
% position of the left camera
[F2, T2, E2, inliers2] = estimateF(pin.c2p, pin.c1, pin.K, 2);  % estimate E/F and decompose the essential
pout.est2.F = F2;
pout.est2.T = T2;
pout.est2.E = E2;
pout.est2.inliers = inliers2;
pout.est2.x1 = pin.c2p;
pout.est2.x2 = pin.c1;
pout.est2.name = 'previous right vs. current left';

T1 = inv(T1);
T2 = inv(T2);

t1 = T1(1:3,4);
R1 = T1(1:3, 1:3);
t2 = T2(1:3,4);
R2 = T2(1:3,1:3);

% Solve for scale of the motion
% A = [t1, -t2];
% b = pin.t0;
% [U,S,V] = svd(A);
% b = U'*b;
% d = diag(S);
% y = nan(size(A,2),1);
% for j=1:length(d)
%     y(j) = b(j)/d(j);
% end
% c = V*y;

P1 = [0 0 0]';
P2 = pin.t0;
c = [-t1'*t1 t2'*t1; -t1'*t2 t2'*t2]\[ P2'*t1-P1'*t1; P2'*t2-P1'*t2];

P1 = P1 - c(1)*t1;
P2 = P2 - c(2)*t2;

t1 = P1 + (P2 - P1)/2;
t2 = t1 - pin.t0;

% update motion direction of the cameras
T1 = [R1 t1; 0 0 0 1];
T2 = [R2 t2; 0 0 0 1];

pout.est1.T_final = T1;
pout.est2.T_final = T2;

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


function [F_est, T_est, E_est, inliers] = estimateF(x1, x2, K, thresh)

if nargin == 3
    thresh = 2; %px
end
% Estimate fundamental: x2'*F*x1
[F_est, inliers] = estimateFundamentalMatrix(x1', x2', 'DistanceType',...
    'Sampson', 'Method', 'RANSAC', 'DistanceThreshold', thresh, 'NumTrials', 500);
inliers = find(inliers);

% compute essential
E_est = K'*F_est*K;

% Decompose essential
T_est = decompose_essential(E_est, K, [x1;x2]);

% refine motion params
x1 = x1(:,inliers);
x2 = x2(:,inliers);
N = size(x1,2);
p1 = nan(3,N);
p2 = nan(3,N);
for i = 1:N
    p1(:,i) = inv(K)*[x1(:,i); 1];
    p2(:,i) = inv(K)*[x2(:,i); 1];
end
[R, t] = essential_lm(T_est(1:3,1:3), T_est(1:3,4), p1, p2);
T_est = [R t; 0 0 0 1];
E_est = skew(t)*R;
F_est = inv(K')*E_est*inv(K);
end

function txt = myupdatefcn(empt, event_obj, x, X)
% Customizes text of data tips

pos = get(event_obj,'Position');
d = x(3:4,:) - repmat(pos',[1 length(x)]);
[~, i] = min(sum(d.*d));
txt = {['Depth: ', num2str(X(3,i))], ['x: ', sprintf('%d %d', x(1,i), x(3,i))], ['d: ', sprintf('x:%d y:%d', x(1,i)-x(3,i), x(2,i)-x(4,i))]};
end
