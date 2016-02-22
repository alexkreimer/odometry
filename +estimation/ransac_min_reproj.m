function pose = ransac_min_reproj(K,base,Xp,xl,xr)

sample_size = 3;
num_pts = length(Xp);
max_inliers = 0;
for i = 1:50
    sample = randsample(num_pts,sample_size)';

    world_points = Xp(:,sample);
    image_points = K\util.e2h(xl(:,sample));
    
	for j = 1:sample_size
        image_points(:,j) = image_points(:,j)/norm(image_points(:,j));
    end
    
    poses = p3p(world_points, image_points);
    max_visible = -1;
    for j = 1:4
        % use chierality to choose the correct solution
        t = poses(1:3,(j-1)*4+1);
        R = poses(1:3,(j-1)*4+2:(j-1)*4+4);
        
        if ~isreal(R) || ~isreal(t)
            continue;
        end

        P1 = K*[eye(3) zeros(3,1)];
        P2 = K*[R' -R'*t];

        [~,visible] = util.triangulate_chieral(xl,xr,P1,P2);
        if sum(visible) > max_visible
            max_visible = sum(visible);
            R0 = R';
            t0 = -R'*t;
        end
    end

    r0 = vrrotmat2vec(R0);
    r0 = r0(1:3)*r0(4);
    rt0 = [r0 t0'];
    
    fun = @(rt) objective(K,base,Xp(:,sample),xl(:,sample),xr(:,sample),rt);
    opt = optimset( optimset('lsqnonlin') ,...
        'Algorithm','levenberg-marquardt',...
        'Diagnostics','off',...
        'Display','off');
    [rt,resnorm,residual,exitflag,output] = lsqnonlin(fun,rt0,[],[],opt);

    residual = objective(K,base,Xp,xl,xr,rt);
    inliers = residual < 1.4*2;
    if sum(inliers) > max_inliers
        best_inliers = inliers;
        best_rt = rt;
    end
end

fun = @(rt) objective(K,base,Xp(:,best_inliers),xl(:,best_inliers),xr(:,best_inliers),rt);
opt = optimset( optimset('lsqnonlin') ,...
    'Algorithm','levenberg-marquardt',...
    'Diagnostics','off',...
    'Display','on',...
    'MaxFunEvals', 100000,...
    'MaxIter', 100000);
[rt,resnorm,residual,exitflag] = lsqnonlin(fun,best_rt,[],[],opt);

% extract parameters
vt = rt(1:3);
R  = vrrotvec2mat(estimation.vtheta2r(vt));
t  = rt(4:6)';

pose = [R t; 0 0 0 1];
end

function res = objective(K,base,Xp,xl,xr,rt)

% extract parameters
vt = rt(1:3);
R  = vrrotvec2mat(estimation.vtheta2r(vt));
t  = rt(4:6)';

% compute reprojection errors
P1 = K*[R t];
P2 = K*[R t-[base;0;0]];

xlp= util.h2e(P1*util.e2h(Xp));
xrp= util.h2e(P2*util.e2h(Xp));

res = sum((xlp-xl).*(xlp-xl)+(xrp-xr).*(xrp-xr),1);

if 0
    % reweight
    weights = exp(-3*abs(res)/max(abs(res)));
    res = res.*weights;
end

end

