function [X,d,J] = triangulate_naive(x1, x2, b, f, cu, cv)
% base is stereo base line
% f is focal distance in pixels
% cu, cv are the coords of the principal point
if iscell(cu)
    cu1 = cu{1};
    cu2 = cu{2};
else
    cu1 = cu;
    cu2 = cu;
end

X = nan(3,1);
J = nan(3,4);

x1cam = x1(1)-cu1;
x2cam = x2(1)-cu2;

%d = max(x1cam(1)-x2cam,0.0001);
d = x1cam-x2cam;
if abs(d) == 0, return; end
X(1) = b*(x1(1)-cu1)/d;
X(2) = b*(x1(2)-cv)/d;
X(3) = b*f/d;

% x = [x1(1) x1(2) x2(1) x2(2)]
% J(i,j) = dX(i)/dx(j)

J(1,1) = b/d-b*x1cam/(d*d); J(1,2) = 0; J(1,3) = b*x1cam/(d*d); J(1,4) = 0;
J(2,1) = -b*(x1(2)-cv)/(d*d); J(2,2) = b/d;  J(2,3) = b*(x1(2)-cv)/(d*d); J(2,4) = 0;
J(3,1) = -b*f/(d*d); J(3,2) = 0; J(3,3) = b*f/(d*d); J(3,4) = 0;

% numeric computation of the Jacobian using derivest
% J0 = jacobianest(@(x) triangulate_naive([x(1) x(2)],[x(3), x(4)],...
% param.base,param.calib.f,param.calib.cu,param.calib.cv), [f1(match(1,j)).pt', f2(match(2,j)).pt']);

end
