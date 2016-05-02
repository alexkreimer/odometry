function [t, predict, inliers] = ransac_minimize_reproj1(X,R,observe,param, varargin)

p = inputParser;
p.addOptional('t0', 2, @isnumeric);
p.addOptional('predict', 0, @isnumeric);
p.KeepUnmatched = true;
parse(p,varargin{:});

% t is our best estiamte
% predict is the predicted location of features under t (computed only for
% inliers)
% inliers are the features that the estimate is computed upon
m = size(observe,1);
status  = false;
num_pts = length(X);
sample_size = 1;
best_inliers = [];

if p.Results.predict
    [rt(1),rt(2),rt(3)] = util.decompose_rotation(R);
    rt(4:6) = p.Results.t0(:);
    [~, residual, predict] = computeJ(X, rt, observe, param, 1:num_pts);
    t = rt(4:6);
    inliers = compute_inliers(residual, param.inlier_thresh, num_pts, m);
    return;
end

for i=1:param.ransac_iter
    sample              = randsample(num_pts,sample_size)';
    [rt(1),rt(2),rt(3)] = util.decompose_rotation(R);
    rt(4:6)             = zeros(1,3);
    for j=1:100
        [J, residual, ~] = computeJ(X,rt,observe,param,sample);
        
        % solving normal equations is a viable way, *but* it has a drawback
        % of possible (unnecessary) conditioning problem, because J'J is
        % computed
        %         JtJ = J'*J;
        %         rc = rcond(JtJ);
        %         if isnan(rc) || rc<1e-12
        %             break;
        %         end
        %         p_gn = (J'*J)\(J'*residual);
        [U,S,V] = svd(J, 'econ');
        c = S(1,1)/S(2,2);
        
        if S(1,1) < 1e-6
            break;
        end
        
        if c>10e4
            warning('J is poorly conditioned');
            break;
        end
        
        p_gn = V*inv(S)*U'*residual;
        if norm(p_gn)<1e-10,
            status = true;
            break;
        else
            rt(4:6) = rt(4:6) + p_gn';
        end
    end
    
    if status,
        [~ ,residual, ~] = computeJ(X,rt,observe,param,1:num_pts);
        inliers = compute_inliers(residual, param.inlier_thresh, num_pts, m);
        if length(inliers) > length(best_inliers)
            best_rt = rt;
            best_inliers = inliers;
        end
    end
end

if isempty(best_inliers)
    return;
end

% final optimization iterations
for k = 1:5
    for i=1:50
        [J, residual, predict] = computeJ(X,best_rt,observe,param,best_inliers);
        
        %     JtJ = J'*J;
        %     rc = rcond(JtJ);
        %     if isnan(rc) || rc<1e-12
        %         break;
        %     end
        %     p_gn = (J'*J)\(J'*residual);
        
        [U,S,V] = svd(J, 'econ');
        c = S(1,1)/S(2,2);
        
        if S(1,1) < 1e-6
            break;
        end
        
        if c>10e4
            warn('J is poorly conditioned');
            break;
        end
        
        p_gn = V*(S\U'*residual);
        
        if norm(p_gn)<1e-10,
            break;
        end
        best_rt(4:6) = best_rt(4:6)+p_gn';
    end
    ninl = length(best_inliers);
    [~, residual, ~] = computeJ(X, best_rt, observe, param, 1:num_pts);
    best_inliers = compute_inliers(residual, param.inlier_thresh, num_pts, m);
    if length(best_inliers)<=ninl
        break;
    end
end

[~, ~, predict] = computeJ(X, best_rt, observe, param, 1:num_pts);
inliers  = best_inliers;
t        = best_rt(4:6)';
end

function [J,residual,predict] = computeJ(X,rt,observe,param,sample)
rx = rt(1); ry = rt(2); rz = rt(3);
tx = rt(4); ty = rt(5); tz = rt(6);

% precompute sine/cosine
sx = sin(rx); cx = cos(rx); sy = sin(ry);
cy = cos(ry); sz = sin(rz); cz = cos(rz);

% compute rotation matrix
r00      = +cy*cz;           r01    = -cy*sz;           r02    = +sy;
r10      = +sx*sy*cz+cx*sz;  r11    = -sx*sy*sz+cx*cz;  r12    = -sx*cy;
r20      = -cx*sy*cz+sx*sz;  r21    = +cx*sy*sz+sx*cz;  r22    = +cx*cy;

m        = size(observe,1);
N        = length(sample);
J        = nan(m*N,3);
predict  = nan(m,N);
residual = nan(m*N,1);

for i = 1:N
    % get 3d point in previous coordinate system
    k   = sample(i);
    X1p = X(1,k);
    Y1p = X(2,k);
    Z1p = X(3,k);
    
    % compute 3d point in current left coordinate system
    X1c = r00*X1p+r01*Y1p+r02*Z1p+tx;
    Y1c = r10*X1p+r11*Y1p+r12*Z1p+ty;
    Z1c = r20*X1p+r21*Y1p+r22*Z1p+tz;
    
    % weighting
    if 1,
        if m == 4
            weight = 1.0/(abs(observe(4*i-3)-param.calib.cu)/abs(param.calib.cu) + 0.05);
        elseif m == 2
            weight = 1.0/(abs(observe(2*i-1)-param.calib.cu)/abs(param.calib.cu) + 0.05);
        end
    else
        weight = 1.0;
    end
    
    % compute 3d point in current right coordinate system
    X2c = X1c-param.base;
    
    % predicted locations for feature point under current [r,t]
    predict(1,i) = param.calib.f*X1c/Z1c+param.calib.cu;
    predict(2,i) = param.calib.f*Y1c/Z1c+param.calib.cv;
    if m == 4
        predict(3,i) = param.calib.f*X2c/Z1c+param.calib.cu;
        predict(4,i) = param.calib.f*Y1c/Z1c+param.calib.cv;
    end
    
    % current residuals
    if m == 2
        residual(2*i-1) = weight*(observe(1,k)-predict(1,i));
        residual(2*i-0) = weight*(observe(2,k)-predict(2,i));
    elseif m == 4
        residual(4*i-3) = weight*(observe(1,k)-predict(1,i));
        residual(4*i-2) = weight*(observe(2,k)-predict(2,i));
        residual(4*i-1) = weight*(observe(3,k)-predict(3,i));
        residual(4*i-0) = weight*(observe(4,k)-predict(4,i));
    end
    
    % iterate over parameters
    for j = 1:3
        switch j
            case 1
                % dr_i/dtx
                X1cd = 1; Y1cd = 0; Z1cd = 0;
            case 2
                % dr_i/dty
                X1cd = 0; Y1cd = 1; Z1cd = 0;
            case 3
                % dr_i/dtz
                X1cd = 0; Y1cd = 0; Z1cd = 1;
        end
        
        % set jacobian entries (project via K)
        if m == 2
            J(2*i-1,j) = weight*param.calib.f*(X1cd*Z1c-X1c*Z1cd)/(Z1c*Z1c);
            J(2*i-0,j) = weight*param.calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c);
        elseif m == 4
            J(4*i-3,j) = weight*param.calib.f*(X1cd*Z1c-X1c*Z1cd)/(Z1c*Z1c);
            J(4*i-2,j) = weight*param.calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c);
            J(4*i-1,j) = weight*param.calib.f*(X1cd*Z1c-X2c*Z1cd)/(Z1c*Z1c);
            J(4*i-0,j) = weight*param.calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c);
        end
    end
end
end

function inliers = compute_inliers(residuals, thresh, num_pts, m)
dist = util.colnorm(reshape(residuals.*residuals, [m num_pts]));
inliers = find(dist < thresh);
end