function [T,F,inliers,residual,success] = rel_motion_F(K,x1,x2)

[F_est,T,~,inliersF,residual,success,~] = estimation.estimateF(x1,x2,K,2);
% WA - we only want negative z here
if T(3,4)>0
    T = inv(T);    
end

if nargout > 1
    inliers = inliersF;
    if nargout > 2
        F = F_est;
    end
end
end
