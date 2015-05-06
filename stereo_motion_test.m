function stereo_motion_test()

dbstop if error;
close all;

param.num = 2;                                 % number of camera motions to generate
param.ad = 6;                                  % length of motion parameterization vector
param.num_pts = 12;                            % number of 3d points to generate
param.K = [718.8560, 0, 607.1928; 0, 718.8560, 185.2157; 0, 0, 1.0000]; % camera intrinsics

% calibration of the stereo rig;
% R0 is the orientation of the right camera in the frame of the left
% t0 is the center of the right camera described in the left frame
R0 = eye(3); t0 = [1,0,0]';

P1 = param.K*[eye(3) zeros(3,1)];              % initial position of the first camera
P2 = param.K*[R0', -R0'*t0];                   % initial position of the 2nd camera in the stereo rig

a = nan(param.ad,param.num);                   % camera parameters
X = 3+rand(3,param.num_pts);                     % 3d points, as seen in the world frame
x = nan(2*param.num,2*param.num_pts);          % image projections, rows 2j,2j+1 hold x,y coordinates for image j
x(1:2,(1:param.num_pts)) = h2e(P1*e2h(X));     % cam1 projection into (left) image
x(1:2,(param.num_pts+1):(2*param.num_pts)) =  h2e(P2*e2h(X)); % cam2 projection into (right) image

a(:,1) = zeros(param.ad,1);                    % camera 1 is located at (0,0,0), its principal axis coincides with z+  direction
% describes translation and orientation of the frame {i} relative to frame {i-1}
a(:,2) = rand(6,1);

T = tr2mat(tget(a,2));

% orientation of frame {2} relative to frame {1}
R = T(1:3,1:3);
    
% origin of frame {2} as seen in frame {1}
t = T(1:3,4);
    
% camera {2} that projects from the world coordinates that coincide with frame {1}
P1_new = param.K*[R',-R'*t];
    
% orientation of the new left frame described in the initial right frame
R2 = R*R0';

% position of the new left origin described in the initial right frame
t2 = R0'*t-R0'*t0;
    
P2_new = param.K*[R2', -R2'*t2];
    
% observations in the current image pair
x(3:4,1:param.num_pts) = h2e(P1_new*e2h(X));
x(3:4,(param.num_pts+1):(2*param.num_pts)) = h2e(P2_new*e2h(X));

% estimate motion
a_est = estimate_stereo_motion(x,param.K,X,param.num_pts,R0,t0);

assert(norm(a_est'-a(:,2),'fro')<1e-12);
end
function a1 = tget(a,j)
% does transformation composition
T = eye(4);

for i=1:j
    T = T*tr2mat(a(:,i));
end

a1 = mat2tr(T);
end

