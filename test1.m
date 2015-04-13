function test1()
close all;
dbstop if error;

KITTI_HOME = '/home/kreimer/tmp/KITTI/dataset';


seq_home = fullfile(KITTI_HOME, 'sequences', '01');
first = 1;
last = 5;

% number of camera motions to generate
param.num = 3;
P1 = nan(3,4,param.num);
P2 = nan(3,4,param.num);
O = nan(param.num,6);

param.lm_max_iter = 200;
param.lm_res_thresh = 1e-5;
param.lm_step_thresh = 1e-5;

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
param.avg_num = 2;
param.num_pts = 3;
param.assert = false;

% All coordinates are described wrt to the Universe Coordinate System or
% some other frame that is described in terms of UCS.

a = repmat({zeros(param.ad,1)},param.num,1);   % camera parameters
X = 20*param.base*rand(param.bd,param.num_pts);% 3d points
x = nan(4,param.num_pts,param.num);            % image projections
x(1:2,:,1) = h2e(P1(:,:,1)*e2h(X));            % cam1 projection into image 1
x(3:4,:,1) = h2e(P2(:,:,1)*e2h(X));            % cam2

% these are used for visualization
O(1,1:3) = [0,0,0];                            % cam1 location relative to {1}
O(1,4:6) = [param.base,0,0];                   % cam2
ax(:,:,1)= eye(3);                             % cam1/2 orientation relative to {1}

% camera 1 is located at (0,0,0), its principal axis coincides with z+
% direction
a{1} = zeros(param.ad,1);
for i=2:param.num
    % orientation of the frame {i} relative to frame {i-1}
    a{i} = [0,0,0,0,0,1.1]';

    T = tr2mat(tget(a,1,i)); R = T(1:3,1:3);
    % origin/axis of the camera {i} relative to {1}
    O(i,1:3) = h2e(T*e2h(O(1,1:3)'))';
    O(i,4:6) = h2e(T*e2h(O(1,4:6)'))';
    ax(:,:,i)= R*ax(:,:,1);
    
    % cameras {i} that projects from frame {1}
    P1(:,:,i) = param.K*[R',-R'*O(i,1:3)'];
    P2(:,:,i) = param.K*[R',-R'*O(i,4:6)'];
    
    % 3d landmarks mapped to frame {i}
    X(:,:,i) = h2e(T\e2h(X(:,:,1)));
    
    % observations in the current image pair
    x(1:2,:,i) = h2e(P1(:,:,i)*e2h(X(:,:,1)));
    x(3:4,:,i) = h2e(P2(:,:,i)*e2h(X(:,:,1)));

    assert(norm(x(1:2,:,i)-h2e(P1(:,:,1)*e2h(X(:,:,i))),'fro')/size(X,2)<1e-6);
    assert(norm(x(3:4,:,i)-h2e(P2(:,:,1)*e2h(X(:,:,i))),'fro')/size(X,2)<1e-6);
    assert(norm(trg(x(:,:,i),param)-X(:,:,i),'fro')/size(X,2)<1e-6);
    A = bsxfun(@minus,X(:,:,i),mean(X(:,:,i),2));
    B = bsxfun(@minus,X(:,:,i-1),mean(X(:,:,i-1),2));
end

figure; plot_axis(O,ax,X);
figure; plot(x(1,:,1),x(2,:,1),'.'); hold on; plot(x(3,:,1),x(4,:,1),'.g'); legend('left view','right view');

sigma = [0,1,2,3,4];
a_est = repmat({nan(6,1)},param.num,length(sigma)); a0_est = a_est;
a_ba = cell(param.num,length(sigma)); a_ba1 = a_ba;
b_ba = cell(param.num_pts,length(sigma)); b_ba1 = b_ba;
a0_ba = cell(size(a_ba));
b0_ba = cell(size(b_ba));
x_noisy = repmat({nan(size(x))},length(sigma),1);
sig0 = repmat({ones(param.obd,1)},param.num_pts,param.num);
for i=1:length(sigma)
    if sigma(i)==0,
        sig1 = sig0; x_noisy{i} = x;
    else
        [sig1,x_noisy{i}] = noise_x(sigma(i),x,param);
    end
    
    a_est{1,i} = a{1}; a0_est{1,i} = a{1};
    for j=2:param.num
        [a_est{j,i},~,a0_est{j,i}] = ransac_minimize_reproj(X(:,:,j-1),x_noisy{i}(:,:,j),param);
        a_est{j,i} = tinv(a_est{j,i});
        a0_est{j,i}= tinv(a0_est{j,i});
    end

    if i==1,
         for j=2:param.num
            % verify that in the noiseless case we are ok
            if param.assert && i==1,
                assert(norm(a_est{j,1}-a{j},'fro')<1e-6);
            end
        end
    end
    
    % BA
    a0_ba{1,i} = a{1};
    for j=2:param.num
        a0_ba{j,i} = tinv(tget(a_est(:,i),1,j));
    end
    
    % convert observations into the cell array
    x_observed = mat2cell(permute(x_noisy{i},[2 3 1]),ones(param.num_pts,1),ones(param.num,1),param.obd);
    x_observed = cellfun(@(x) permute(x,[3 2 1]), x_observed, 'UniformOutput', false);
    
    % initial guess for the structure
    b0_ba(:,i) = mat2cell(permute(X(:,:,1),[3 2 1]),1,ones(param.num_pts,1),3)';
    b0_ba(:,i) = cellfun(@(X) permute(X,[3 2 1]), b0_ba(:,i), 'UniformOutput', false);
    
    if param.assert && i == 1,
        for j=2:param.num
            T = tr2mat(a0_ba{j,i});
            Xt = h2e(T*e2h(X(:,:,1)));
            %assert(norm(Xt-X(:,:,j),'fro')<1e-6);
            xt = h2e(P1(:,:,1)*e2h(Xt));
            assert(norm(xt-x(1:2,:,j),'fro')/size(x,2)<1e-3);
            xt = h2e(P2(:,:,1)*T*e2h(X(:,:,1)));
            assert(norm(xt-x(3:4,:,j),'fro')/size(x,2)<1e-3);
        end
    end
    
    % run bundle adjustment with dummy covariance
    [a_ba(:,i),b_ba(:,i)] = lmlsq(x_observed,sig0,a0_ba(:,1),b0_ba(:,1),@obj_ba,param);

    % run BA with real covariance
    [a_ba1(:,i),b_ba1(:,i)] = lmlsq(x_observed,sig1,a0_ba(:,1),b0_ba(:,1),@obj_ba,param);

end
    
dt = zeros(8,length(sigma));
for i=1:length(sigma)
    for j=1:param.num
        [d1,d2] = tdist(tinv(tget(a,1,j)), tinv(tget(a_est(:,i),1,j)));
        dt(1,i) = dt(1,i)+d1;
        dt(2,i) = dt(2,i)+d2;
        [d1,d2] = tdist(tinv(tget(a,1,j)), a_ba{j,i});
        dt(3,i) = dt(3,i)+d1;
        dt(4,i) = dt(4,i)+d2;
        [d1,d2] = tdist(tinv(tget(a,1,j)), a_ba1{j,i});        
        dt(5,i) = dt(5,i)+d1;
        dt(6,i) = dt(6,i)+d2;
        [d1,d2] = tdist(tinv(tget(a,1,j)), tinv(tget(a0_est(:,i),1,j)));
        dt(7,i) = dt(7,i)+d1;
        dt(8,i) = dt(8,i)+d2;
    end
end
dt = dt/param.num;
h = figure();
set(h,'defaulttextinterpreter','latex');
subplot(221); plot(sigma,dt(1,:),sigma,dt(7,:)); legend('Subsequent frame reprojection minimization','Orthogonal procrustes'); xlabel('$\sigma_{max}$'); ylabel('$\left\|\theta_{ML}-\theta\right\|_2^2$'); title('Subsequent frame parameter estimation');
subplot(222); plot(sigma,dt(3,:),sigma,dt(5,:)); legend('BA rot error', 'BAwC rot error');  xlabel('$\sigma_{max}$'); ylabel('$\left\|\theta_{ML}-\theta\right\|_2^2$'); title('BA parameter refinement');
subplot(223); plot(sigma,dt(2,:),sigma,dt(8,:)); legend('Subsequent frame reprojection minimization','Orthogonal procrustes'); xlabel('$\sigma_{max}$'); ylabel('$\left\|\theta_{ML}-\theta\right\|_2^2$');  title('Subsequent frame parameter estimation');
subplot(224); plot(sigma,dt(4,:),sigma,dt(6,:)); legend('BA translation error','BAwC translation error');  xlabel('$\sigma_{max}$'); ylabel('$\left\|\theta_{ML}-\theta\right\|_2^2$'); title('BA parameter refinement');
end

function [sig,x_noisy] = noise_x(max_std,x,param)
    sig = repmat({zeros(param.obd)},param.num_pts,param.num);
    x_noisy = nan(size(x));
    for i=1:size(x,1)
        for j=1:size(x,2)
            for k=1:size(x,3)
                sig{j,k} = randi(max_std*max_std,param.obd,1);
                R = chol(diag(sig{j,k}));
                x_noisy(:,j,k) = x(:,j,k) + R*randn(param.obd,1);
            end
        end
    end
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
plot3(X(1,:,1),X(2,:,1),X(3,:,1),'.');
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
end

function [dr,dt] = tdist(a1,a2)
T1 = tr2mat(a1);
T2 = tr2mat(a2);
R1 = T1(1:3,1:3);
R2 = T2(1:3,1:3);
t1 = T1(1:3,4);
t2 = T2(1:3,4);

dr = norm(a1(1:3)-a2(1:3)); dt = norm(t1-t2);
end



function a1 = tget(a,i1,i2)
% this function assumes that a{j} holds transformation parameters that
% describe the pose of frame {j} relative to frame {j-1}

T = eye(4);
for i=i1:i2
    T = T*tr2mat(a{i});
end
a1 = mat2tr(T);
end