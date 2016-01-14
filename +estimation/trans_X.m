function T = trans_X(K,R,base,Xp,x1,x2)

fun = @(t) objective(K,R,base,Xp,x1,x2,t);
[t_opt,resnorm,residual,exitflag] = lsqnonlin(fun,zeros(3,1));
T = [R t_opt; 0 0 0 1];
end

function val = objective(K,R,base,Xp,x1,x2,t)

val = nan(1,length(Xp));

for i=1:length(Xp)
    P1 = K*[R t];
    P2 = K*[R t+R'*[base 0 0]'];
    
    x1p= util.h2e(P1*util.e2h(Xp(:,i)));
    x2p= util.h2e(P2*util.e2h(Xp(:,i)));
    
    val(i) = norm(x1p-x1(:,i))+norm(x2p-x2(:,i));
end
end
