function sparse_lba(x,sigma,tr,X,visible,param)

% each 2 columns of x correspond to view i (left/right)
% each 2 rows of x correspond to feature j
% the values should be nan if feature j is not visible in view j

% column i of tr should contain the initial pose guess for view i
% column i of X should contain the initial guess for 3d point X_i


end

function ba_J(x,sigma,tr,visible,X,param)
% frames.visible - index into X of features visible in this frame
% frames.x - observations, 4 x visible_num (1:2 left frame, 3:4 right frame)
% frames.tr - frame pose, 6 x 1 (3 for euler angles, 3 for translation)
% frames.A - predictor derivative w.r.t. motion params, 4*size(frame.x,1) x 6
% frames.B - predictor derivative w.r.t. structure params, 4*size(frame.x,1) x 3*visible_num
% X 3d points

% number of 3d points
num_pts = size(X,2);
% number of views
num_frm = length(frames);
A = cell(num_pts,num_frm);
predict = nan(num_pts,num_frm,4);
r = nan(num_pts,num,frm,4);

for j = 1:num_frm
    % motion parameteres, rotation(euler angles)/translation
    rx = frames(j).tr(1); ry = frames(j).tr(2); frames(j).rz = tr(3);
    tx = frames(j).tr(4); ty = frames(j).tr(5); frames(j).tz = tr(6);
    
    % precompute sine/cosine
    sx = sin(rx); cx = cos(rx); sy = sin(ry); cy = cos(ry); sz = sin(rz); cz = cos(rz);
    
    % compute rotation matrix
    r00    = +cy*cz;           r01    = -cy*sz;           r02    = +sy;
    r10    = +sx*sy*cz+cx*sz;  r11    = -sx*sy*sz+cx*cz;  r12    = -sx*cy;
    r20    = -cx*sy*cz+sx*sz;  r21    = +cx*sy*sz+sx*cz;  r22    = +cx*cy;
    
    % the following 2 loops iterate over all visible features across all
    % images
    for i=1:num_pts
        % preallocate
        A{i,j} = nan(4,6);
        
        % if point i is not visible in view j: A_ij = 0
        if visible{i,j} == false
            A{i,j} = zeros(4,6);
            continue;
        end
        
        X1p = X(1,i); Y1p = X(2,i); Z1p = X(3,i);
        
        % [R t]*e2h(X)
        X1c = r00*X1p+r01*Y1p+r02*Z1p+tx;
        Y1c = r10*X1p+r11*Y1p+r12*Z1p+ty;
        Z1c = r20*X1p+r21*Y1p+r22*Z1p+tz;
        
        % predictions h2e(K[R t]e2h(X))
        predict(i,j,1) = param.calib.f*X1c/Z1c+param.calib.cu;
        predict(i,j,2) = param.calib.f*Y1c/Z1c+param.calib.cv;
        predict(i,j,3) = param.calib.f*X2c/Z1c+param.calib.cu;
        predict(i,j,4) = param.calib.f*Y1c/Z1c+param.calib.cv;
        
        r(i,j,1) = x(i,j,1) - predict(i,j,1);
        r(i,j,2) = x(i,j,2) - predict(i,j,2);
        r(i,j,3) = x(i,j,3) - predict(i,j,3);
        r(i,j,4) = x(i,j,3) - predict(i,j,4);
        
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
        
        % derivation w.r.t. motion parameters
        for p=1:6
            switch p
                case 1
                    % dR/drx*X
                    X1cd = 0;
                    Y1cd = rdrx10*X1p+rdrx11*Y1p+rdrx12*Z1p;
                    Z1cd = rdrx20*X1p+rdrx21*Y1p+rdrx22*Z1p;
                case 2
                    % dR/dry*X
                    X1cd = rdry00*X1p+rdry01*Y1p+rdry02*Z1p;
                    Y1cd = rdry10*X1p+rdry11*Y1p+rdry12*Z1p;
                    Z1cd = rdry20*X1p+rdry21*Y1p+rdry22*Z1p;
                case 3
                    % dR/drz*X
                    X1cd = rdrz00*X1p+rdrz01*Y1p;
                    Y1cd = rdrz10*X1p+rdrz11*Y1p;
                    Z1cd = rdrz20*X1p+rdrz21*Y1p;
                case 4
                    % dR/dtx*X
                    X1cd = 1; Y1cd = 0; Z1cd = 0;
                case 5
                    % dR/dty*X
                    X1cd = 0; Y1cd = 1; Z1cd = 0;
                case 6
                    % dR/dtz*X
                    X1cd = 0; Y1cd = 0; Z1cd = 1;
            end
            % division by Z1c*Z1c is because of the derivative of h2e
            A{i,j}(1,p) = param.calib.f*(X1cd*Z1c-X1c*Z1cd)/(Z1c*Z1c);
            A{i,j}(2,p) = param.calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c);
            A{i,j}(3,p) = param.calib.f*(X1cd*Z1c-X2c*Z1cd)/(Z1c*Z1c);
            A{i,j}(4,p) = param.calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c);
        end
    end
end


B = cell(num_pts,num_frm);
% iterate over all structure parameters
for i=1:num_pts
    for j=1:num_frm
        % motion parameteres, rotation(euler angles)/translation
        rx = frames(j).tr(1); ry = frames(j).tr(2); frames(j).rz = tr(3);
        tx = frames(j).tr(4); ty = frames(j).tr(5); frames(j).tz = tr(6);
        
        % precompute sine/cosine
        sx = sin(rx); cx = cos(rx); sy = sin(ry);
        cy = cos(ry); sz = sin(rz); cz = cos(rz);
        
        % compute rotation matrix
        r00 = +cy*cz;          r01 = -cy*sz;          r02 = +sy;
        r10 = +sx*sy*cz+cx*sz; r11 = -sx*sy*sz+cx*cz; r12 = -sx*cy;
        r20 = -cx*sy*cz+sx*sz; r21 = +cx*sy*sz+sx*cz; r22 = +cx*cy;
        
        % number of visible 3d points should be the same as the number of
        % observed features provided
        Xv = X(:,frames(j).visible);
        X1p = Xv(1,j); Y1p = Xv(2,j); Z1p = Xv(3,j);
        % [R t]*e2h(X)
        X1c = r00*X1p+r01*Y1p+r02*Z1p+tx;
        Y1c = r10*X1p+r11*Y1p+r12*Z1p+ty;
        Z1c = r20*X1p+r21*Y1p+r22*Z1p+tz;
        
        % derivation w.r.t structure parameters (3d points)
        for p=1:3
            switch j
                case 1
                    X1cd = r00; Y1cd = r10; Z1cd = r20;
                case 2
                    X1cd = r01; Y1cd = r11; Z1cd = r21;
                case 3
                    X1cd = r02; Y1cd = r12; Z1cd = r22;
            end
            B{i,j}(1,p) = param.calib.f*(X1cd*Z1c-X1c*Z1cd)/(Z1c*Z1c);
            B{i,j}(2,p) = param.calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c);
            B{i,j}(3,p) = param.calib.f*(X1cd*Z1c-X2c*Z1cd)/(Z1c*Z1c);
            B{i,j}(4,p) = param.calib.f*(Y1cd*Z1c-Y1c*Z1cd)/(Z1c*Z1c);
        end
    end
end

% precompute U blocks
U = cell(length(frames),1);
for j=1:length(frames)
    U{j} = zeros(6);
    for i=1:size(X,2)
        if visible{i,j}
            U{j} = U{j} + A{i,j}'*sigma{i,j}*A{i,j};
        end
    end
end

% precompute V blocks
V = cell(size(X,2),1);
for i=1:size(X,2)
    for j=1:length(frames)
        if visible{i,j},
            V{j} = V{j} + B{i,j}'*sigma{i,j}*B{i,j};
        end
    end
end

% precompute W
W = cell(length(frames),size(X,2));
for i=1:length(frames)
    for j=1:size(X,2)
        if visible{i,j}, W{i,j} = A{i,j}'*sigma{i,j}*B{i,j};
        else W{i,j} = zeros(6,3); end
    end
end

ra = cell(length(frames),1);
for j = 1:length(frames)
    for i=1:size(X,2)
        if visible{i,j},
            ra{j} = ra{j} + A{i,j}'*sigma{i,j}*r{i,j};
        end
    end
end

rb = cell(size(X,2),1);
for i=1:size(X,2)
    for j=1:length(frames)
        if visible{i,j}
            rb{i} = rb{i} + B{i,j}'*sigma{i,j}*r{i,j};
        end
    end
end
end
