function [R,inliers,residual] = rel_motion_H(K,x1,x2,d,b,varargin)

% x1 interest points in image 1
% x2 interest points in image 2
% d  depth of the points
% b  baseline
% size(x1)=size(x2)=[2,N]
p = inputParser;
p.KeepUnmatched = true;

p.addOptional('depth_thr',100,@isnumeric);
parse(p,varargin{:});
depthThr = p.Results.depth_thr;

assert(all(size(x1)==size(x2)));
[rows, ~] = size(x1);
if rows==2
    x1 = util.e2h(x1);
    x2 = util.e2h(x2);
end
inliers = d>depthThr*b;
if sum(inliers)<10
    warning('not enough distant points for current distance threshold')
    inliers = d>.7*depthThr*b;
    if sum(inliers)<10
        inliers = d>.5*depthThr*b;
    end
end

x1o = x1(:,~inliers);
x2o = x2(:,~inliers);
x1  = x1(:, inliers);
x2  = x2(:, inliers);
[H,R,inliers_H,residual] = estimation.H_inf_nonlin(K,x1,x2,varargin{:});
%t = estimation.trans_geom(K,H,x1o,x2o);
t = [0;0;0];
T = [R t; 0 0 0 1];
inliers = find(inliers);
inliers = inliers(inliers_H);
end