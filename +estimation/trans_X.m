function t = trans_X(K,R,base,Xp,xl,xr,t0)

if nargin < 7
    t0 = zeros(3,1);
end

sample_size = 1;
num_pts = length(Xp);
max_inliers = 0;
for i = 1:500
    sample = randsample(num_pts,sample_size);
    fun = @(t) objective(K,R,base,Xp(:,sample),xl(:,sample),xr(:,sample),t);
    options = optimset( optimset('lsqnonlin'), 'Algorithm','levenberg-marquardt','Display','off');
    [t_opt,resnorm,residual,exitflag] = lsqnonlin(fun,t0,[],[],options);
    
    residual = objective(K,R,base,Xp,xl,xr,t_opt);

    r1 = residual(1:num_pts);
    r2 = residual((num_pts+1):end);

    inliers = r1<1.4*1 & r2<1.4*1;
    n = sum(inliers);
    if n > max_inliers
        best_inliers = inliers;
        best_t = t_opt;
        if n==num_pts
            break;
        end
    end
end

fun = @(t) objective(K,R,base,Xp(:,best_inliers),xl(:,best_inliers),xr(:,best_inliers),t);
options = optimset( optimset('lsqnonlin') ,...
    'Algorithm','levenberg-marquardt',...
    'Display','on',...
    'MaxFunEvals', 10000,...
    'MaxIter', 10000);
[t,resnorm,residual,exitflag] = lsqnonlin(fun,best_t,[],[],options);
end

function res = objective(K,R,base,Xp,xl,xr,t)
P1 = K*[R t];
P2 = K*[R t-[base;0;0]];

xlp= util.h2e(P1*util.e2h(Xp));
xrp= util.h2e(P2*util.e2h(Xp));

res = sum(sqrt([(xlp-xl).*(xlp-xl) (xrp-xr).*(xrp-xr)]),1);

if 0
    % reweight
    weights = exp(-3*abs(res)/max(abs(res)));
    res = res.*weights;
end

end
