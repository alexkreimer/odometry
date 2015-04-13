function [a,b] = lmlsq(x,sigma,a,b,predict,param)
% x{i,j} are the coordinates of point i in view j. if the cell x{i,j} is
% empty it means that point i is not visible in view j.
% sigma{i,j} is the covariance of x{i,j}
% a{j} are camera parameters for view j
% b{i} are 3d point i coordinates
% the procedure returns the refined parameters

s = cell2mat(sigma');
s = inv(diag(s(:)));

lambda = 1e-2; ni = 10; status = false;
for i=1:param.lm_max_iter
    [r,J] = predict(x,a,b,param);
    
    if r'*s*r < param.lm_res_thresh,
        fprintf('iter %d obj rms:%g, norm(Jr)=%g\n',i,sqrt(r'*r/length(r)),norm(J'*r));
        break;
    end
    while 1,
        JtJ = J'*s*J;
        % JtJ = JtJ + lambda*diag(diag(JtJ));
        JtJ = JtJ + lambda*eye(size(JtJ));
        if rcond(JtJ)<1e-15 || isnan(rcond(JtJ)),
            lambda = lambda*ni;
            continue;
        end
        delta = -inv(JtJ)*(J'*s)*r;
        if norm(delta)<param.lm_step_thresh,
            status = true;
            fprintf('done after %d iterations %g\n',i);
            break;
        end
        a_upd = cell2mat(a);
        a_upd = a_upd + delta(1:size(a_upd,1));
        a_upd = mat2cell(a_upd,param.ad*ones(size(a,1),1));
        
        b_upd = cell2mat(b);
        b_upd = b_upd + delta(param.ad*size(a,1)+1:end);
        b_upd = mat2cell(b_upd,param.bd*ones(size(b,1),1));
        [r_upd, ~] = predict(x,a_upd,b_upd,param);

        if r'*r > r_upd'*r_upd
            a = a_upd;
            b = b_upd;
            lambda = lambda/ni;
            break;
        else
            lambda = lambda*ni;
        end
    end
    fprintf('iter %d, obj rms:%g,',i,sqrt(r'*r/length(r)));
    fprintf('delta: %g, lambda=%g, norm(Jr)=%g\n',norm(delta),lambda,norm(J'*r));
    if status, break, end;
end

end

