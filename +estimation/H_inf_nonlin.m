function [H,R_out,inliers_out,residual_out] = H_inf_nonlin(K,x1,x2,varargin)

% x1 distant interest points in image 1
% x2 distant interest points in image 2
% d  depth of the points
% b  baseline
% size(x1)=size(x2)=[2,N]
p = inputParser;

p.addOptional('absRotInit',false);
p.addOptional('F',[]);
p.addOptional('inlier_thr',1,@isnumeric);
p.addOptional('ransac_iter',100,@isnumeric);

p.KeepUnmatched = true;

parse(p,varargin{:});

assert(all(size(x1)==size(x2)));

[rows, ~] = size(x1);

if rows==2
    x1 = util.e2h(x1);
    x2 = util.e2h(x2);
end

% Estimate fundamental: x2'*F*x1. Normalized 8-point algorithm
num_pts  = length(x1);
s        = 3;
max_size = -1;
thr2     = p.Results.inlier_thr;

opt = optimset(optimset('lsqnonlin'),...
    'Algorithm','levenberg-marquardt',...
    'Diagnostics','off','Display','off');

for i = 1:p.Results.ransac_iter
    sample = randsample(num_pts,s);
    p1 = normc(K\x1(:,sample));
    p2 = normc(K\x2(:,sample));
    
    if p.Results.absRotInit
        [regParams,~,~] = util.absor(p1,p2,'doTrans',false);
        r = vrrotmat2vec(regParams.R);
        vt0 = r(1:3)*r(4);
    else
        vt0 = zeros(1,3);
    end

%     if ~isempty(p.Results.F)
%         fun = @(vt) objective1(K,p.Results.F,x1(:,sample),x2(:,sample),vt);
%     else
%         fun = @(vt) objective2(K,x1(:,sample),x2(:,sample),vt);
%     end
% 
%     [vt,resnorm,result,exitflag,output] = lsqnonlin(fun,vt0,[],[],opt);
%     R = vrrotvec2mat(estimation.vtheta2r(vt));
    %fprintf('homography symm reprojection error RMS is %g\n',sqrt(d*d'/length(d)));
%     H = K*R/K;
    H = K*regParams.R/K;
    inliers = util.homogdist2d(H,x1,x2,thr2);
    support_size = length(inliers);
    %fprintf('support size: %d\n',support_size);
    if support_size>max_size
        max_size = support_size;
        inliers_best = inliers;
        %vt_best = vt;
        vt_best = vt0;
    end
end

if ~isempty(p.Results.F)
    fun = @(vt) objective1(K,p.Results.F,x1(:,inliers_best),x2(:,inliers_best),vt);
else
    fun = @(vt) objective2(K,x1(:,inliers_best),x2(:,inliers_best),vt);
end
opt = optimset(optimset('lsqnonlin') , 'Algorithm','levenberg-marquardt', 'Diagnostics','off', 'Display','off');
[vt,resnorm,residual,exitflag] = lsqnonlin(fun,vt_best,[],[],opt);
R = vrrotvec2mat(estimation.vtheta2r(vt));
H = K*R/K;
%fprintf('homography symm reprojection error RMS is %g\n',sqrt(d*d'/length(d)));

% Hx1 = util.h2e(H*x1);
% plot([x2(1,:); Hx1(1,:)],[x2(2,:); Hx1(2,:)],'-om', 'DisplayName', 'Refined homography');
% legend show;
if nargout>1
    R_out = R;
    if nargout > 2
        inliers_out = inliers;
        if nargout>3
            residual_out = residual;
        end
    end
end
end

function val = objective1(K,F,x1,x2,vtheta)
r = estimation.vtheta2r(vtheta);
R = vrrotvec2mat(r);
H = K*R/K;

e1 = util.h2e(x2)-util.h2e(H*x1);
e2 = util.h2e(x1)-util.h2e(H\x2);
if norm(F)>1e-5
    l1 = F*x1;
    l1 = l1./repmat(util.colnorm(l1(1:2,:)),[3 1]);
    e1_ortho = diag(l1(1:2,:)'*e1);
    
    
    l2 = F'*x2;
    l2 = l2./repmat(util.colnorm(l2(1:2,:)),[3 1]);
    e2_ortho = diag(l2(1:2,:)'*e2);
    
    val = [e1_ortho.*e1_ortho; e2_ortho.*e2_ortho]';
else
    val = [util.colnorm(e1); util.colnorm(e2)]';
end
end

function val = objective2(K,x1,x2,vtheta)

% convert v*\theta into [v \theta]
vec = estimation.vtheta2r(vtheta);
R = vrrotvec2mat(vec);

% H_inf
H = K*R/K;

% residuals
e1 = util.h2e(x2)-util.h2e(H*x1);
e2 = util.h2e(x1)-util.h2e(H\x2);

% objective value
val = [util.colnorm(e1); util.colnorm(e2)]';
end

function val = objective3(K,F,x1,x2,vtheta)
r = estimation.vtheta2r(vtheta);
R = vrrotvec2mat(r);
H = K*R/K;

ex1 = util.h2e(x1);
ex2 = util.h2e(x2);

Hx1 = util.h2e(H*x1);
Hix2= util.h2e(H\x2);

figure; hold on;
plot([ex1(1,:); Hix2(1,:)],[ex1(2,:); Hix2(2,:)],'b--o');

e1 = util.h2e(x2)-util.h2e(H*x1);
e2 = util.h2e(x1)-util.h2e(H\x2);

epipole  = util.h2e(null(F));
plot(epipole(1),epipole(2),'r*');

epilines = epipolarLine(F,ex1');
points = lineToBorderPoints(epilines, [376 1241]);
line(points(:, [1,3])', points(:, [2,4])');

epipole  = util.h2e(null(F'));
plot(epipole(1),epipole(2),'r*');

epilines = epipolarLine(F',ex1');
points = lineToBorderPoints(epilines, [376 1241]);
line(points(:, [1,3])', points(:, [2,4])');

if norm(F)>1e-5
    l1 = F*x1;
    l1 = l1./repmat(util.colnorm(l1(1:2,:)),[3 1]);
    e1_ortho = diag(l1(1:2,:)'*e1);

    l2 = F'*x2;
    l2 = l2./repmat(util.colnorm(l2(1:2,:)),[3 1]);
    e2_ortho = diag(l2(1:2,:)'*e2);
    
    val = [e1_ortho.*e1_ortho; e2_ortho.*e2_ortho]';
else
    val = [util.colnorm(e1); util.colnorm(e2)]';
end
end

