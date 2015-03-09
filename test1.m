function test1()
close all;
dbstop if error;

KITTI_HOME = '/home/kreimer/tmp/KITTI/dataset';


seq_home = fullfile(KITTI_HOME, 'sequences', '01');
first = 1;
last = 5;

% number of camera motions to generate
param.num = 5;
P1 = nan(3,4,param.num);
P2 = nan(3,4,param.num);
O = nan(param.num,6);

% setup camera parameters (KITTI)
P1(:,:,1) = [718.8560, 0, 607.1928,0; 0, 718.8560, 185.2157,0; 0, 0, 1.0000, 0];
P2(:,:,1) = [718.8560, 0, 607.1928, -386.1448; 0, 718.8560, 185.2157, 0; 0, 0, 1.0000, 0];

% baseline, focal, principal point
param.base = -P2(1,4,1)/P1(1,1,1);
param.calib.f = P1(1,1,1); 
param.calib.cu = P1(1,3,1);
param.calib.cv = P1(2,3,1);
param.K = P1(1:3,1:3,1);
% this commands minimization procedure to first solve 3d-3d rigid motion
% and use it as an initial point for Gauss-Newton reprojection minimization 
% procedure
param.init = true;
% fundamental
param.F = vgg_F_from_P(P1(:,:,1),P2(:,:,1));
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
num_pts = 10;

% camera 1 is located at (0,0,0) its principal axis coincides with z+
% direction

% camera motions
a = repmat({zeros(6,1)},param.num,1);

X = 20*param.base*rand(3,num_pts);
x = nan(4,num_pts,param.num);

x(1:2,:,1) = h2e(P1(:,:,1)*e2h(X)); x(3:4,:,1) = h2e(P2(:,:,1)*e2h(X));

% origin and axis of the first pose
O(1,1:3) = [0,0,0]; O(1,4:6) = [param.base,0,0]; ax(:,:,1) = eye(3);

for i=2:param.num
    % camera motion
    a{i} = [0,0,0,0,0,1.1]'; T = tr2mat(a{i}); R = T(1:3,1:3);
    % origin/axis of the 2nd pose
    O(i,1:3) = h2e(T*e2h(O(i-1,1:3)'))'; O(i,4:6) = h2e(T*e2h(O(i-1,4:6)'))'; ax(:,:,i) = R*ax(:,:,1);
    % camera matrices of the 2nd pair
    P1(:,:,i) = param.K*[R',-R'*O(i,1:3)']; P2(:,:,i) = param.K*[R',-R'*O(i,4:6)'];
    % observations in the 2nd image pair
    x(1:2,:,i) = h2e(P1(:,:,i)*e2h(X)); x(3:4,:,i) = h2e(P2(:,:,i)*e2h(X));
end

figure; plot_axis(O,ax,X);
figure; plot(x(1,:,1),x(2,:,1),'.'); hold on; plot(x(3,:,1),x(4,:,1),'.g'); legend('left view','right view');

sigma = linspace(0,4,100);
d = nan(2,length(sigma));
a_est = repmat({nan(6,1)},param.num,length(sigma));
for i=1:length(sigma)
    for j=2:param.num
        N = 10;
        for k=1:N
            a_est{j,i}(:,k) = ransac_minimize_reproj(X,x(:,:,j)+sigma(i)*randn(size(x(:,:,j))),param);
        end
    end
end

figure; plot(sigma,d(1,:),sigma,d(2,:)); legend('rotation error', 'translation error');
end

function plot_axis(O,ax,X)
% plot axis/structure
grid on; hold on; axis on; daspect([1,1,1]); view([-37 30]);
for i=1:size(O,1)
    plot3(O(i,1),O(i,2),O(i,3),'or');  plot3(O(i,4),O(i,5),O(i,6),'or');
    arrow3(repmat(O(i,1:3),3,1),repmat(O(i,1:3),3,1)+ax(:,:,i),'b');
    arrow3(repmat(O(i,4:6),3,1),repmat(O(i,4:6),3,1)+ax(:,:,i),'g');
    text(O(i,1)+ax(1,1,i),O(i,2)+ax(1,2,i),O(i,3)+ax(1,3,i),'x');
    text(O(i,1)+ax(2,1,i),O(i,2)+ax(2,2,i),O(i,3)+ax(2,3,i),'y');
    text(O(i,1)+ax(3,1,i),O(i,2)+ax(3,2,i),O(i,3)+ax(3,3,i),'z');
end
plot3(X(1,:),X(2,:),X(3,:),'.');
axis equal;
end
