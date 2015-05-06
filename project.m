function x = project(X,K,pose)
% transform into the camera frame
X = pose*[X;ones(1,size(X,2))];
x = K*h2e(X);
x(3,x(3,:)<0) = nan;
end