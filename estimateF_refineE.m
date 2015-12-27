function [F, T, E, inliers, F_lin] = estimateF_refineE(x1, x2, K, thresh)

if nargin == 3
    thresh = 2; %px
end

% Estimate fundamental: x2'*F*x1. Normalized 8-point algorithm
[F, inliers] = estimateFundamentalMatrix(x1', x2',...
    'DistanceType', 'Sampson',...
    'Method', 'RANSAC',...
    'DistanceThreshold',...
    thresh, 'NumTrials', 500);
inliers = find(inliers);

F_lin = F;

% refine motion params
x1 = e2h(x1(:,inliers));
x2 = e2h(x2(:,inliers));



end

