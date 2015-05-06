function alpha = armijos(FUN,x,d)

[fval, grad] = FUN(x);
max_iter = 100;
sigma = 1e-4;
alpha = 4;
beta  = .02;
c = grad'*d;
k = 1;

% alpha(1) = alpha;
% for j=2:100
%     alpha(j) = alpha(j-1)*beta;
%     [fval(j),~] = FUN(x+alpha(j)*d);
% end
% 
% figure; plot(alpha,fval); title('phi(alpha)');
best_delta = inf;

while k < max_iter
    [fval1, ~] = FUN(x+alpha*d);
    delta = fval1-fval;
    if delta<best_delta
        best_alpha = alpha;
    end
    
    if fval1 - fval < alpha*sigma*c;
        return;
    else
        alpha = beta*alpha;
    end
    k = k + 1;
end

warning('armijo failed');
alpha = 0;

end