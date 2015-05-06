function s = invE(E,P,K,sol)
%CentralCamera.invE Decompose essential matrix
%
% S = C.invE(E) decomposes the essential matrix E (3x3) into the camera motion.
% In practice there are multiple solutions and S (4x4xN) is a set of homogeneous
% transformations representing possible camera motion.
%
% S = C.invE(E, P) as above but only solutions in which the world point P is visible
% are returned.
%
% Reference::
% Hartley & Zisserman,
% "Multiview Geometry",
% Chap 9, p. 259
%
% Y.Ma, J.Kosecka, S.Soatto, S.Sastry,
% "An invitation to 3D",
% Springer, 2003.
% p116, p120-122
%
% Notes::
% - The transformation is from view 1 to view 2.
%
% See also CentralCamera.E.

% we return T from view 1 to view 2

[U,S,V] = svd(E);
% singular values are (sigma, sigma, 0)
if nargin<4
    sol = 0;
end

if sol
    % H&Z solution
    W = [0 -1 0; 1 0 0; 0 0 1];   % rotz(pi/2)
    
    t = U(:,3);
    R1 = U*W*V';
    if det(R1) < 0,
        disp('flip');
        V = -V;
        R1 = U*W*V';
        det(R1)
    end
    R2 = U*W'*V';
    
    % we need to invert the solutions since our definition of pose is
    % from initial camera to the final camera
    s(:,:,1) = inv([R1 t; 0 0 0 1]);
    s(:,:,2) = inv([R1 -t; 0 0 0 1]);
    s(:,:,3) = inv([R2 t; 0 0 0 1]);
    s(:,:,4) = inv([R2 -t; 0 0 0 1]);
    p1 = project(P,K,s(:,:,1));
    p2 = project(P,K,s(:,:,2));
    p3 = project(P,K,s(:,:,3));
    p4 = project(P,K,s(:,:,4));
else
    % Ma etal solution, p116, p120-122
    % Fig 5.2 (p113), is wrong, (R,t) is from camera 2 to 1
    if det(V) < 0
        V = -V;
        S = -S;
    end
    if det(U) < 0
        U = -U;
        S = -S;
    end
    R1 = U*rotz(90)'*V';
    R2 = U*rotz(-90)'*V';
    t1 = vex(U*rotz(90)*S*U');
    t2 = vex(U*rotz(-90)*S*U');
    % invert (R,t) so its from camera 1 to 2
    s(:,:,1) = inv( [R1 t1; 0 0 0 1] );
    s(:,:,2) = inv( [R2 t2; 0 0 0 1] );
end

if nargin > 2
    for i=1:size(s,3)
        if ~any(isnan(project(P,K,s(:,:,i))))
            s = s(:,:,i);
            fprintf('solution %d is good\n', i);
            return;
        end
    end
    warning('no solution has given point in front of camera');
end
end