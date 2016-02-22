function [t,residual,inliers] = ransac_minimize_reproj1(X,R,xl,xr,param)

status  = false;
num_pts = length(X);
sample_size = 1;
best_inliers = [];
for i=1:500
    sample       = randsample(num_pts,sample_size)';
    [rt(1),rt(2),rt(3)] = util.decompose_rotation(R);
    rt(4:6)             = zeros(1,3);
    for j=1:100
        [J, residual, ~] = computeJ(X,rt,[xl;xr],param,sample);
        JtJ = J'*J;
        rc = rcond(JtJ);
        if isnan(rc) || rc<1e-12
            break;
        end
        p_gn = (J'*J)\(J'*residual);
        if norm(p_gn)<1e-10,
            status = true;
            break;
        else
            rt(4:6) = rt(4:6) + p_gn';
        end
    end
    
    if status,
        sample = 1:size(X,2);
        [~ ,residual, ~] = computeJ(X,rt,[xl;xr],param,sample);
        inliers = compute_inliers(residual,4,sample);
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
for i=1:50
    [J, residual, predict] = computeJ(X,best_rt,[xl;xr],param,best_inliers);
    JtJ = J'*J;
    rc = rcond(JtJ);
    if isnan(rc) || rc<1e-12
        break;
    end
    p_gn = (J'*J)\(J'*residual);
    if norm(p_gn)<1e-10,
        break;
    end
    best_rt(4:6) = best_rt(4:6)+p_gn';
end
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
r00    = +cy*cz;           r01    = -cy*sz;           r02    = +sy;
r10    = +sx*sy*cz+cx*sz;  r11    = -sx*sy*sz+cx*cz;  r12    = -sx*cy;
r20    = -cx*sy*cz+sx*sz;  r21    = +cx*sy*sz+sx*cz;  r22    = +cx*cy;

J = nan(4*length(sample),3);
predict = nan(4,length(sample));
residual = nan(4*length(sample),1);

for i = 1:length(sample)
    % get 3d point in previous coordinate system
    X1p = X(1,sample(i));
    Y1p = X(2,sample(i));
    Z1p = X(3,sample(i));
    
    % compute 3d point in current left coordinate system
    X1c = r00*X1p+r01*Y1p+r02*Z1p+tx;
    Y1c = r10*X1p+r11*Y1p+r12*Z1p+ty;
    Z1c = r20*X1p+r21*Y1p+r22*Z1p+tz;
    
    % weighting
    if 1,
        weight = 1.0/(abs(observe(4*i-3)-param.calib.cu)/abs(param.calib.cu) + 0.05);
    else
        weight = 1.0;
    end
    
    % compute 3d point in current right coordinate system
    X2c = X1c-param.base;
    
    % predicted locations for feature point under current [r,t]
    predict(1,i) = param.calib.f*X1c/Z1c+param.calib.cu; % left u
    predict(2,i) = param.calib.f*Y1c/Z1c+param.calib.cv; % left v
    predict(3,i) = param.calib.f*X2c/Z1c+param.calib.cu; % right u
    predict(4,i) = param.calib.f*Y1c/Z1c+param.calib.cv; % right v
    
    % current residuals
    residual(4*i-3) = weight*(observe(1,sample(i))-predict(1,i));
    residual(4*i-2) = weight*(observe(2,sample(i))-predict(2,i));
    residual(4*i-1) = weight*(observe(3,sample(i))-predict(3,i));
    residual(4*i-0) = weight*(observe(4,sample(i))-predict(4,i));
    
    for j=1:3
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
        J(4*i-3,j) = weight*param.calib.f*(X1cd*Z1c-X1c*Z1cd)/(Z1c*Z1c); % left u'
        J(4*i-2,j) = weight*param.calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c); % left v'
        J(4*i-1,j) = weight*param.calib.f*(X1cd*Z1c-X2c*Z1cd)/(Z1c*Z1c); % right u'
        J(4*i-0,j) = weight*param.calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c); % right v'
    end
end
end

function inliers = compute_inliers(residuals, thresh, active)
dist = reshape(residuals.*residuals, [4 numel(active)]);
inliers = dist < thresh*thresh;
inliers = inliers(1,:) & inliers(2,:);
inliers = find(inliers);
end