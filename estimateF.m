function [F, T, E, inliers,F_lin] = estimateF(x1, x2, K, thresh, ref)

if nargin < 5
    ref = 'tangent';
end
if nargin == 3
    thresh = 2; %px
end

% Estimate fundamental: x2'*F*x1. Normalized 8-point algorithm
[F, inliers] = estimateFundamentalMatrix(x1', x2',...
    'DistanceType', 'Sampson',...
    'Method', 'RANSAC',...
    'DistanceThreshold',...
    thresh, 'NumTrials', 500);
inliers = find(inliers);

% output
F_lin = F; E = nan(3); T = nan(4);

if strcmp(ref, 'none')
    return;
end

% refine motion params
x1 = e2h(x1(:,inliers));
x2 = e2h(x2(:,inliers));
options = optimset( optimset('lsqnonlin') ,...
    'Algorithm','levenberg-marquardt',...
    'Diagnostics','off',...
    'TolFun',1e-8,...
    'TolX',1e-10,...
    'Display','off');
    %'MaxFunEvals', 20000,...
    %MaxIter', 10000);

if strcmp(ref,'direct') || strcmp(ref,'all')
    fun = @(f) sampsonF(reshape(f,[3 3]), x1, x2);
    [f, resnorm,residual,exitflag] = lsqnonlin(fun, F_lin(:), [], [], options);

    if exitflag < 1
        disp('direct LM refinement of the fundamental failed');
    elseif strcmp(ref,'all')
        F(:,:,1) = reshape(f,[3 3]);
        E(:,:,1) = K'*F*K;
        T(:,:,1) = [decompose_essential(E, K, [x1;x2]); 0 0 0 1];
    else
        F = reshape(f,[3 3]);
        E = K'*F*K;
        T = [decompose_essential(E, K, [x1;x2]); 0 0 0 1];        
    end
end

if strcmp(ref,'vtheta') || strcmp(ref,'all')
    c0  = vtheta_init_params(F_lin,K,x1,x2);
    fun = @(c) vtheta_sampson(c, K, x1, x2);
    [c, resnorm, residual,exitflag] = lsqnonlin(fun, c0, [], [], options);

    if exitflag < 1
        disp('vtheta LM refinement of the essential failed');
    elseif strcmp(ref,'all')
        [F(:,:,2),E(:,:,2),T(:,:,2)] = vtheta_param2F(c,K);
    else
        [F,E,T] = vtheta_param2F(c,K);        
    end
end

if strcmp(ref, 'tangent') || strcmp(ref, 'all')
    T_temp = decompose_essential(K'*F_lin*K, K, [x1;x2]);
    R = T_temp(1:3,1:3);
    t = T_temp(1:3,4);

    h0  = quaternion.rotationmatrix(R).e;
    c0  = [zeros(3,1); t];
    fun = @(c) tangent_sampson(c(1:6), h0, K, x1, x2);
    [c, resnorm, residual,exitflag] = lsqnonlin(fun, c0, [], [], options);

    if exitflag < 1
        disp('tangent LM refinement of the essential failed');
    elseif strcmp(ref,'all')
        [F(:,:,3),E(:,:,3),T(:,:,3)] = tangent_param2F(c,h0,K);
    else
        [F,E,T] = tangent_param2F(c,h0,K);
    end
end

end

% tangent plane parameterization support functions
function [F,E,T] = tangent_param2F(c,h0,K)
v = c(1:3);
t = c(4:6);
hZ = tangent_param2quaternion(v, h0);
R = quaternion(hZ).RotationMatrix;
E = skew(t)*R;
F = K'\E/K;
T = [R t; 0 0 0 1];
end

function val = tangent_sampson(c,h0,K,x1,x2)
% extract rotation and translation parameters
[F,~,~] = tangent_param2F(c,h0,K);
val = sampsonF(F, x1, x2);
end

function hZ = tangent_param2quaternion(v, h0)
% convert 3-vector v and the operating point h0 into a quaternion hZ

B = tangent_ortho_basis(h0);

% convert parameter vector to 4-vector and normalize it
v4  = B*v;

if norm(v4) == 0
    hZ = h0;
else
    v4N = v4/norm(v4);
    
    % resulting quaternion
    theta = norm(v4);
    hZ = sin(theta)*v4N + cos(theta)*h0;
end

end

% produce orthogonal basis to the radius one sphere in R4
function B = tangent_ortho_basis(h)
% note, norm(h)=1 so one of the ifs has to happen
if h(1) ~= 0
    b1 = [-h(2)/h(1); 1; 0; 0];
    b2 = [-h(3)/h(1); 0; 1; 0];
    b3 = [-h(4)/h(1); 0; 0; 1];
elseif h(2) ~= 0
    b1 = [1; -h(1)/h(2); 0; 0];
    b2 = [0; -h(3)/h(2); 1; 0];
    b3 = [0; -h(4)/h(2); 0; 1];
elseif h(3) ~= 0
    b1 = [1; 0; -h(1)/h(3); 0];
    b2 = [0; 1; -h(2)/h(3); 0];
    b3 = [0; 0; -h(4)/h(3); 1];
elseif h(4) ~= 0
    b1 = [1; 0; 0; -h(1)/h(4)];
    b2 = [0; 1; 0; -h(2)/h(4)];
    b3 = [0; 0; 1; -h(3)/h(4)];
end
B       = [b1 b2 b3];
[U,~,~] = svd(B);
B       = U(:,1:3);
end

% v*theta rotation parameterization support functions
function [F,E,T] = vtheta_param2F(c,K)
    angle = norm(c(1:3));
    if angle == 0
        axis = [0;0;1];
    else
        axis  = c(1:3)/angle;
    end
    q  = quaternion.angleaxis(angle,axis);
    R  = q.RotationMatrix;
    t = c(4:6);
    T = [R t; 0 0 0 1];
    E = skew(t)*R;
    F = K'\E/K;
end
    
function val = vtheta_sampson(c,K,x1,x2)
    [F,~,~] = vtheta_param2F(c,K);
    val = sampsonF(F,x1,x2);
end

function c = vtheta_init_params(F,K,x1,x2)
    E = K'*F*K;
    T = decompose_essential(E, K, [x1;x2]);
    R = T(1:3,1:3);
    t = T(1:3,4);
    q = quaternion.rotationmatrix(R);
    [axis, angle] = q.AngleAxis();
    c = [axis*angle; t];
end
