function val = objective1(w,K,base,tr1,tr2,ratio,sigma,c)

err1 = objective1_sampson(tr1,base,K,c(1));
err2 = objective1_sampson(tr2,base,K,c(2));

delta = (abs(ratio)-(1+abs(c(1)/c(2))));

n1 = length(err1);
n2 = length(err2);

val = [err1/n1; err2/n2; w*delta/(sigma*sigma)];
end

function err = objective1_sampson(tr,base,K,c)

% linear solution inliers
x1 = tr.x1;
x2 = tr.x2;

T = inv(tr.T);
R = T(1:3,1:3);
t = c(1)*(T(1:3,4)+[base;0;0])-[base;0;0];
F = K'\(util.skew(t)*R)/K;
err = estimation.sampsonF(F, util.e2h(x1), util.e2h(x2));
end
