function [best_tr, best_inliers, tr0, predict, rms] = ransac_minimize_reproj(X, observe, param)
%% RANSAC & LM to estimate rigid motion parameters
% _X_ is 3d point cloud, s.t. size(X) = [3 N]
%
% observe are projections of X in a new frame (i.e., after camera
% rotation/translation).  The function assumes that
% observed(1:2,j)/observe(3:4,j) are the projections of X(:,j).  It also
% allows for some degree of mismatches.

% param is a structure that holds various parameters (TODO: describe
% params here)
% best_tr is a 6-vector of motion params that will be estimated.
% best_inliers is a support set for best_tr (i.e., best_tr is estimated
% based on best_inliers)
% tr0 is an initial guess that was used for LM (it is obtained by solving
% 3d-2-3d rigid motion)
% predict(:,j) is a predicted projection of point X(:,j) into the current
% image plane
% rms is a reprojection error rms

status = false;
best_tr = zeros(6,1);
tr0 = [];
predict = [];
rms = 0;
best_inliers = [];
tr0 = nan(6,1);
for j=1:param.ransac_iter
    active = datasample(1:size(X,2),param.model_size,'Replace',false);
    tr = zeros(6,1);
    if param.init
        X0 = X(:,active); X1 = trg(observe(:,active),param);
        if ~any(isnan(X1(:))) && ~any(isinf(X1(:))) && ~any(isinf(X0(:)))
            tr = solve3d23d(X0,X1); tr0 = tr;
        end
    else
        tr0 = nan(6,1);
    end

    for i=1:param.lm_max_iter
        [J, residual, ~] = computeJ(X, tr, observe, param, active);
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
            tr = tr+p_gn;
        end
    end
    
    if status,
        active = 1:size(X,2);
        [~ ,residual, ~] = computeJ(X, tr, observe, param, active);
        inliers = compute_inliers(residual, param.inlier_thresh, active);
        if length(inliers) > length(best_inliers)
            best_tr = tr;
            best_inliers = inliers;
        end
    end
end

if isempty(best_inliers)
    return;
end

% final optimization iterations
for i=1:50
    [J, residual, predict] = computeJ(X, best_tr, observe, param, best_inliers);
    JtJ = J'*J;
    rc = rcond(JtJ);
    if isnan(rc) || rc<1e-12
        break;
    end
    p_gn = (J'*J)\(J'*residual);
    if norm(p_gn)<1e-10,
        break;
    end
    best_tr = best_tr+p_gn;
end

inliers = compute_inliers(residual, param.inlier_thresh, best_inliers);

predict = reshape(predict, [4 numel(best_inliers)]);
rms = sqrt(sum(residual.*residual)/(4*numel(best_inliers)));

err = colnorm(predict-observe(:,best_inliers));
% histogram(err);
% title('reprojection error distribution');
end

function [J,residual,predict] = computeJ(X,tr,observe,param,active)

% motion parameteres, rotation/translation
rx = tr(1); ry = tr(2); rz = tr(3);
tx = tr(4); ty = tr(5); tz = tr(6);

% precompute sine/cosine
sx = sin(rx); cx = cos(rx); sy = sin(ry);
cy = cos(ry); sz = sin(rz); cz = cos(rz);

% compute rotation matrix
r00    = +cy*cz;           r01    = -cy*sz;           r02    = +sy;
r10    = +sx*sy*cz+cx*sz;  r11    = -sx*sy*sz+cx*cz;  r12    = -sx*cy;
r20    = -cx*sy*cz+sx*sz;  r21    = +cx*sy*sz+sx*cz;  r22    = +cx*cy;
% dR/drx
rdrx10 = +cx*sy*cz-sx*sz;  rdrx11 = -cx*sy*sz-sx*cz;  rdrx12 = -cx*cy;
rdrx20 = +sx*sy*cz+cx*sz;  rdrx21 = -sx*sy*sz+cx*cz;  rdrx22 = -sx*cy;
% dR/dry
rdry00 = -sy*cz;           rdry01 = +sy*sz;           rdry02 = +cy;
rdry10 = +sx*cy*cz;        rdry11 = -sx*cy*sz;        rdry12 = +sx*sy;
rdry20 = -cx*cy*cz;        rdry21 = +cx*cy*sz;        rdry22 = -cx*sy;
% dR/drz
rdrz00 = -cy*sz;           rdrz01 = -cy*cz;
rdrz10 = -sx*sy*sz+cx*cz;  rdrz11 = -sx*sy*cz-cx*sz;
rdrz20 = +cx*sy*sz+sx*cz;  rdrz21 = +cx*sy*cz-sx*sz;

J = nan(4*length(active),6);
predict = nan(4,length(active));
residual = nan(4*length(active),1);

for i = 1:length(active)
    % get 3d point in previous coordinate system
    X1p = X(1,active(i));
    Y1p = X(2,active(i));
    Z1p = X(3,active(i));
    
    % fprintf('observation %d; X1p,X2p,Z1p=%g,%g,%g\n',i,X1p,Y1p,Z1p);
    % compute 3d point in current left coordinate system
    X1c = r00*X1p+r01*Y1p+r02*Z1p+tx;
    Y1c = r10*X1p+r11*Y1p+r12*Z1p+ty;
    Z1c = r20*X1p+r21*Y1p+r22*Z1p+tz;
    % fprintf('observation  X1c,Y1c,Z1c=%g,%g,%g\n',X1c,Y1c,Z1c);
    
    % weighting
    if 1,
        weight = 1.0/(abs(observe(4*i-3)-param.calib.cu)/abs(param.calib.cu) + 0.05);
    else
        weight = 1.0;
    end
    
    % compute 3d point in current right coordinate system
    X2c = X1c-param.base;
    
    % predicted locations for current feature point under current [r,t]
    predict(1,i) = param.calib.f*X1c/Z1c+param.calib.cu; % left u
    predict(2,i) = param.calib.f*Y1c/Z1c+param.calib.cv; % left v
    predict(3,i) = param.calib.f*X2c/Z1c+param.calib.cu; % right u
    predict(4,i) = param.calib.f*Y1c/Z1c+param.calib.cv; % right v
    
    % current residuals
    residual(4*i-3) = weight*(observe(1,active(i))-predict(1,i));
    residual(4*i-2) = weight*(observe(2,active(i))-predict(2,i));
    residual(4*i-1) = weight*(observe(3,active(i))-predict(3,i));
    residual(4*i-0) = weight*(observe(4,active(i))-predict(4,i));
    
    % fprintf('X2c: %g\n',X2c);
    % r[i] = observed-K[R(r),t]*Xp
    % for all paramters do
    for j=1:6
        switch j
            case 1
                % dr_i/drx
                X1cd = 0;
                Y1cd = rdrx10*X1p+rdrx11*Y1p+rdrx12*Z1p;
                Z1cd = rdrx20*X1p+rdrx21*Y1p+rdrx22*Z1p;
            case 2
                % dr_i/dry
                X1cd = rdry00*X1p+rdry01*Y1p+rdry02*Z1p;
                Y1cd = rdry10*X1p+rdry11*Y1p+rdry12*Z1p;
                Z1cd = rdry20*X1p+rdry21*Y1p+rdry22*Z1p;
            case 3
                % dr_i/drz
                X1cd = rdrz00*X1p+rdrz01*Y1p;
                Y1cd = rdrz10*X1p+rdrz11*Y1p;
                Z1cd = rdrz20*X1p+rdrz21*Y1p;
            case 4
                % dr_i/dtx
                X1cd = 1; Y1cd = 0; Z1cd = 0;
            case 5
                % dr_i/dty
                X1cd = 0; Y1cd = 1; Z1cd = 0;
            case 6
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

function tr = solve3d23d(X0,X1)
tr = nan(6,1);
B = bsxfun(@minus,X1,mean(X1,2))';
A = bsxfun(@minus,X0,mean(X0,2))';
[~, R, ~] = procrust(A,B);
[tr(1),tr(2),tr(3)] = decompose_rotation(R);
tr(4:6) = -R*mean(X0,2)+mean(X1,2);
end

function inliers = compute_inliers(residuals, thresh, active)

dist = reshape(residuals.*residuals, [4 numel(active)]);

inliers = dist < thresh*thresh;
inliers = inliers(1,:) & inliers(2,:);
inliers = find(inliers);
end