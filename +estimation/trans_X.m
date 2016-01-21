function t = trans_X(K,R,base,Xp,xl,xr,t0)

if nargin < 7
    t0 = zeros(3,1);
end

fun = @(t) objective(K,R,base,Xp,xl,xr,t);
options = optimset( optimset('lsqnonlin') ,...
    'Algorithm','levenberg-marquardt',...
    'Diagnostics','off',...
    'TolFun',1e-8,...
    'TolX',1e-10,...
    'Display','off',...
    'MaxFunEvals', 3000);
    %MaxIter', 10000);
[t_opt,resnorm,residual,exitflag] = lsqnonlin(fun,t0,[],[],options);
t = t_opt/norm(t_opt);
end

function res = objective(K,R,base,Xp,xl,xr,t)

R = R';
t = -R*t;

P1 = K*[R' -R'*t];
P2 = K*[R' -R'*(t+[base;0;0])];

xlp= util.h2e(P1*util.e2h(Xp));
xrp= util.h2e(P2*util.e2h(Xp));

res = sum((xlp-xl).*(xlp-xl)+(xrp-xr).*(xrp-xr),1);

if 0
% reweight
weights = exp(-3*abs(res)/max(abs(res)));
res = res.*weights;
end

end
