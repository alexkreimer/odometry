function test_epipole_sensitivity(num_pts, z_max, z_min)

OUT_DIR = '/home/kreimer/prj/odometry/debug';

num_frames = 2;
K = [718.8560,  0,      607.1928;
       0,     718.8560, 185.2157;
       0,       0,        1.0000];       % camera intrinsics

% R0/t0 is the orientation/origin of the right camera as seen from the left camera
baseline            = [.5;0;0];
Id                  = eye(3);
T0                  = [Id, baseline; 0 0 0 1];
param.base          = norm(baseline);
param.calib.f       = K(1);
param.calib.cu      = K(1,3);
param.calib.cv      = K(2,3);
param.ransac_iter   = 1000;
param.model_size    = 3;
param.init          = true;
param.lm_max_iter   = 1000;
param.inlier_thresh = 3;
P1                  = K*[Id zeros(3,1)];     % left camera

X(1,:) = -3 + (3-(-3))*rand(1, num_pts); % 3d points, as seen in the world frame
X(2,:) = -3 + (3-(-3))*rand(1, num_pts);
X(3,:) = z_min + (z_max-z_min)*rand(1, num_pts);

x = nan(2*num_frames, 2*num_pts);     % image projections

T = [Id zeros(3,1); 0 0 0 1];
dt = [Id [0;0;1]; 0 0 0 1];
for frame=1:num_frames
    % observations in the new image pair
    
    [x((2*frame-1):2*frame, 1:num_pts), visible1]              = project(P1, X, T);
    [x((2*frame-1):2*frame,(num_pts+1):(2*num_pts)), visible2] = project(P1, X, T*T0);
    
    T = dt*T;
    assert(all([visible1,visible2]),'some points are invisible in the new pose of a rig');
end

% number of experiments
N = 10;

% preallocate
e = nan(2,N,2);


figure;
hold on;
title('feature tracks (both original and noisy)');

for i=1:N
    fprintf('experiment %d\n', i);
    x1 = x(1:2, 1:num_pts);
    x2 = x(3:4, 1:num_pts);
    
    plot([x1(1,:); x2(1,:)], [x1(2,:); x2(2,:)]);
    
    if i>1
        x1 = x1 + randn(size(x1));
        x2 = x2 + randn(size(x2));
    end
    
    plot([x1(1,:); x2(1,:)], [x1(2,:); x2(2,:)]);
    

    
    [F, T, E, inliers] = estimateF(x2, x1, K, 2);
    [U,S,V] = svd(F);
    e(:, i, 1) = h2e(V(:,3));
    err(:,i) = sampson_error(E, K\e2h(x1), K\e2h(x2));
    
    e(:,i,2) = epipole_from_H(x1, x2);    
end

%save_dbg(fullfile(OUT_DIR, sprintf('tracks_%d_%d_%d', num_pts, z_max, z_min)));
%close;

edist = @(x1,x2) min([1, abs(x1-x2)/min([abs(x1),abs(x2)])]);
for i=1:N
    for j = 1:2
        x1 = e(1,1,j); x2 = e(1,i,j);
        y1 = e(2,1,j); y2 = e(2,i,j);
        
        de(1,i,j) = edist(x1, x2);
        de(2,i,j) = edist(y1, y2);
    end
end
err = colnorm(err);

% relative error plot
figure;
hold on;
title('epipole relative coordinate error');
plot(de(:,:,1)');
legend('x','y');

%save_dbg(fullfile(OUT_DIR, sprintf('e_rel_error_%d_%d_%d', num_pts, z_max, z_min)));
%close;

% scatter the epipoles
figure; hold on;
plot(e(1,:,1), e(2,:,1), '.r');
plot(e(1,:,2), e(2,:,2), '.b');
plot(K(1,3),K(2,3), '*g');
legend('epipole(s) from F', 'epipole(s) from H', 'real epipole');

%save_dbg(fullfile(OUT_DIR, sprintf('epipoles_%d_%d_%d', num_pts, z_max, z_min)));
%close;

% fit normal distribution to the epipole estimations
data1 = e(:,:,1)';
data2 = e(:,:,2)';

pd1 = fitdist(data1(:,1), 'Normal');
pd2 = fitdist(data1(:,2), 'Normal');
pd3 = fitdist(data2(:,1), 'Normal');
pd4 = fitdist(data2(:,2), 'Normal');

x_values1 = min(data1(:,1)):.1:max(data1(:,1));
x_values2 = min(data1(:,2)):.1:max(data1(:,2));
x_values3 = min(data2(:,1)):.1:max(data2(:,1));
x_values4 = min(data2(:,2)):.1:max(data2(:,2));

y1 = pdf(pd1, x_values1);
y2 = pdf(pd2, x_values2);
y3 = pdf(pd3, x_values3);
y4 = pdf(pd4, x_values4);

figure; hold on; title('epipole per-coordinate distributions');
plot(x_values1, y1,'LineWidth',2);
plot(x_values2, y2,'LineWidth',2);
plot(x_values3, y3,'LineWidth',2);
plot(x_values4, y4,'LineWidth',2);
legend('x of epiple from F', 'y of epipole from F', 'x of epipole from H', 'y of epipole from H');

%save_dbg(fullfile(OUT_DIR, sprintf('e_pdf_%d_%d_%d', num_pts, z_max, z_min)));
%close;

fileID = fopen(fullfile(OUT_DIR, sprintf('e_pdf_%d_%d_%d.txt', num_pts, z_max, z_min)),'w');

fprintf(fileID, 'x-distribution-F mean: %g, sigma: %g\n', pd1.mu, pd1.sigma);
fprintf(fileID, 'y-distribution-F mean: %g, sigma: %g\n', pd2.mu, pd2.sigma);
fprintf(fileID, 'x-distribution-H mean: %g, sigma: %g\n', pd3.mu, pd3.sigma);
fprintf(fileID, 'y-distribution-H mean: %g, sigma: %g\n', pd4.mu, pd4.sigma);

fclose(fileID);
end

function val = sampson_error(E, x1, x2)
% x1, x2 inlier coordinates, s.t. x2'Ex1=0

N = size(x1,2);

val = nan(N,1);
for i=1:N
    % current point
    p1 = x1(:,i);
    p2 = x2(:,i);
    
    if length(p1) == 2
        p1 = [p1;1];
    end
    
    if length(p2) == 2
        p2 = [p2;1];
    end

    p2tE = p2'*E;
    Ep1  = E*p1;
    
    val(i) = p2'*E*p1/sqrt(p2tE(1)^2 + (p2tE(1))^2 + (Ep1(1))^2 + (Ep1(2))^2);
end
end