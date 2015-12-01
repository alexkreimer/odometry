function test_objective1()

% number of steps to produce data for
steps = 3;

% stereo rig internals/externals
baseline = 0.53;
T0 = [eye(3) [baseline 0 0]'; 0 0 0 1];
T0inv = inv(T0);

K = [707.0912,   0,      601.8873;
       0,      707.0912  183.1104;
       0,        0,        1.0000];

% preallocate
npts = 1000;
P0 = nan(3,4,steps);
P1 = nan(3,4,steps);
x_left  = nan(2,npts,steps);
x_right = nan(2,npts,steps);

%initial cameras   
P0(:,:,1) = K*[eye(3) zeros(3,1)];
P1(:,:,1) = K*T0inv(1:3, :);

% test points
X = rand(3,npts);
X(3,:) = X(3,:) + 10;

% constant camera motion
R = eye(3); t = [0 0 1]'; delta = [R t; 0 0 0 1];

% generate data; initial camera is located at the origin of the world
% coordinates

for i=1:steps
    if i>1
        T(:,:,i) = delta*T(:,:,i-1);
    else
        T(:,:,1) = [eye(3) zeros(3,1); 0 0 0 1];
    end
    
    Tleft  = inv(T(:,:,i));
    Tright = inv(T0*T(:,:,i));
    
    P0(:,:,i) = K*Tleft(1:3,:);
    P1(:,:,i) = K*Tright(1:3,:);
    
    % project
    x_left(:,:,i)  = h2e(P0(:,:,i)*e2h(X));
    x_right(:,:,i) = h2e(P1(:,:,i)*e2h(X));
    
    % add noise
    x_left(:,:,i)  = x_left(:,:,i)  + randn(size(x_left(:,:,i)));
    x_right(:,:,i) = x_right(:,:,i) + randn(size(x_right(:,:,i)));
end

% first 2 motions
T1 = delta;
T2 = delta;

x1 = [K\e2h(x_left(:,:,2)); K\e2h(x_right(:,:,1))];
x2 = [K\e2h(x_left(:,:,3)); K\e2h(x_right(:,:,2))];

ratio = 2+.1*randn; sigma = 1;
fun = @(c) objective1(T0(1:3,4), T1, T2, x1, x2, ratio, sigma, c(1), c(2));

[c, fval] = fminsearch(fun, [1, 1]);
end

