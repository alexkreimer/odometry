function [x,visible] = project(P,X,pose)

if nargin>2
    if size(pose,1) == 3
        pose = [pose;0 0 0 1];
    end
    
    % convert the points to new camera frame given by 'pose'
    X = pose\util.e2h(X);
else
    X = util.e2h(X);
end
x = P*X;

% visibility of the 3d points
if nargout>1
    visible = X(4,:).*x(3,:)>0;
end

x = util.h2e(x);

end