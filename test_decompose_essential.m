function test_decompose_essential()

dbstop if error; close all;

% number of 3d points to generate
num_pts = 12;

% camera intrinsics
% K = [718.8560, 0, 607.1928; 0, 718.8560, 185.2157; 0, 0, 1.0000]; 

K = [645.24,0,635.96;0,645.24,194.13;0,0,1]; % KITTI demo.cpp

% initial camera
P1 = K*[eye(3) zeros(3,1)]; 

X = 20+10*rand(3,num_pts);              % 3d points, as seen in the world frame

x = nan(4,num_pts);                     % image projections
[x(1:2,:),visible] = project(P1,X);     % cam1 projection

assert(all(visible),'all points must be visible in view 1');

% describes translation and orientation of the frame {i} relative to frame {i-1}
a = rand(6,1);
%a = [pi/2,0,0,0,0,1]';

T = tr2mat(a);

% orientation of frame {2} relative to frame {1}
R = T(1:3,1:3);
    
% origin of frame {2} as seen in frame {1}
t = T(1:3,4);
    
% camera {2} that projects from the world coordinates that coincide with frame {1}
P2 = K*[R',-R'*t];
    
% observations in the current image pair
[x(3:4,:),visible] = project(P2,X);
assert(all(visible),'all points must be visible in view 2');

% Estimate 1st essential: x2'*F*x1
x1 = e2h(x(1:2,:));
x2 = e2h(x(3:4,:));
E = K'*fundmatrix(x1,x2)*K; 
% Decompose essential
T1 = decompose_essential(E,K,x);

assertT(T1,T,'inv')
fprintf('test is successfull\n');
end
