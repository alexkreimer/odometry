function gen_test()

dbstop if error;
close all;

param.num = 2;                                 % number of camera motions to generate
param.ad = 6;                                  % length of motion parameterization vector
param.num_pts = 12;                            % number of 3d points to generate
param.K = [718.8560, 0, 607.1928; 0, 718.8560, 185.2157; 0, 0, 1.0000]; % camera intrinsics

% calibration of the stereo rig;
% R0 is the orientation of the right camera in the frame of the left
% t0 is the center of the right camera described in the left frame
R0 = eye(3); t0 = [1,0,0]';

P1 = param.K*[eye(3) zeros(3,1)];              % initial position of the first camera
P2 = param.K*[R0', -R0'*t0];                   % initial position of the 2nd camera in the stereo rig

a = nan(param.ad,param.num);                   % camera parameters
X = 3+rand(3,param.num_pts);                     % 3d points, as seen in the world frame
x = nan(2*param.num,2*param.num_pts);          % image projections, rows 2j,2j+1 hold x,y coordinates for image j
x(1:2,(1:param.num_pts)) = h2e(P1*e2h(X));     % cam1 projection into (left) image
x(1:2,(param.num_pts+1):(2*param.num_pts)) =  h2e(P2*e2h(X)); % cam2 projection into (right) image

a(:,1) = zeros(param.ad,1);                    % camera 1 is located at (0,0,0), its principal axis coincides with z+  direction
for i=2:param.num
    % describes translation and orientation of the frame {i} relative to frame {i-1}
    a(:,i) = [0,.5,0,0,0,1]';
    a(:,i) = rand(6,1);
    
    T = tr2mat(tget(a,i));
    
    % orientation of frame {i} relative to frame {1}
    R = T(1:3,1:3);
    
    % origin of frame {i} as seen in frame {1}
    t = T(1:3,4);
    
    % camera {i} that projects from frame {1}
    P1_new = param.K*[R',-R'*t];
    
    % orientation of the new left frame described in the initial right frame
    R2 = R*R0';
    % position of the new left origin described in the initial right frame
    t2 = R0'*t-R0'*t0;
    
    P2_new = param.K*[R2', -R2'*t2];
    
    % observations in the current image pair
    x((2*i-1):(2*i),1:param.num_pts) = h2e(P1_new*e2h(X));
    x((2*i-1):(2*i),(param.num_pts+1):(2*param.num_pts)) = h2e(P2_new*e2h(X));
    
    a_est = estimate_motion(x,param.K,X,param.num_pts,R0,t0);
    assert(norm(a_est'-a(:,2),'fro')<1e-12);
    %     a0 = rand(6,1);
    %     aopt = a(:,i);
    %     a = newton(x,param.K,a0,aopt);
    
    a0 = a_est';
    aopt = a(:,i);
    a = gauss_newton(x,param.K,a0,aopt);
    
    %     a0 = a0/norm(a0(4:6))+[0 0 0 0 0 .1]';
    %     a = minsearch(x,param.K,a0);
    1;
end
end

function a1 = tget(a,j)
% does transformation composition
T = eye(4);

for i=1:j
    T = T*tr2mat(a(:,i));
end

a1 = mat2tr(T);
end

function a = newton(x,K,a,aopt)

x1 = h2e(K\[x(1:2,:);ones(1,size(x,2))]);
x2 = h2e(K\[x(3:4,:);ones(1,size(x,2))]);
x = [x1;x2];

N = 10;
stat = nan(2,N);
lambda = 1;

% alpha = linspace(1,norm(a-aopt),100);
% u = aopt-a; u = u/norm(u);
% for j=1:100
%     [f,grad,H] = symess(x,a+alpha(j)*u,lambda);
%     v(1,j) = f;
%     v(2,j) = acos(grad'*u/norm(grad))*180/pi;
%     pk = -H\grad;
%     v(3,j) = acos(pk'*u/norm(pk))*180/pi;
% end
%
% figure; plot(v(1,:)); legend('value');
% figure;
% plot(v(2,:)); hold on; plot(v(3,:));
% legend('steepest descent vs. u','newton vs. u');

for j=1:N
    [f,grad,H] = sym_tr_obj(x,a,lambda);
    if norm(grad) < 1e-8
        break;
    end
    % GN search direction
    pk = -H\grad;
    
    % backtracking step size search
    alpha = armijos(@(y) sym_tr_obj(x,y,lambda),a,pk);
    
    if alpha == 0
        break;
    else
        % update
        a = a+pk;
    end
    stat(1,j) = f;
    stat(2,j) = norm(grad);
end
figure; plot(stat(1,:)); title('objective values');
figure; plot(stat(2,:)); title('gradient norms');

end

function a = gauss_newton(x,K,a,aopt)

x1 = h2e(K\[x(1:2,:);ones(1,size(x,2))]);
x2 = h2e(K\[x(3:4,:);ones(1,size(x,2))]);
x = [x1;x2];

N = 1000;        % max iterations
stat = nan(2,N); % stats
lambda = 10;     % weight of the t't-1 constraint
lm = .01;        % initial LM damping parameter
beta = 10;       % LM step parameter

alpha = linspace(1,norm(a-aopt),100);
u = aopt-a;
u = u/norm(u);
for j=1:100
    [r,J] = gaussnewton_tr_obj(x,a+alpha(j)*u,lambda);
    v(1,j) = r'*r;
    g = J'*r;
    v(2,j) = acos(g'*u/norm(g))*180/pi;
    H = J'*J;
    pk = -H\(J'*r);
    v(3,j) = acos(pk'*u/norm(pk))*180/pi;
end

figure; plot(v(1,:)); legend('value');
figure;
plot(v(2,:)); hold on; plot(v(3,:));
legend('grad vs. u','GN vs. u');

k=1;
for j=1:N
    fprintf('\nLM iter %d:',j);
    [r,J] = gaussnewton_tr_obj(x,a,lambda);
    fprintf('obj: %g, gnorm: %g,',r'*r,norm(J'*r));
    if norm(J'*r) < 1e-8
        break;
    end
    
    H = J'*J + lm*eye(6);
    % GN search direction
    pk = -H\(J'*r);
    %pk = u;
    
    % backtracking step size search
    alpha = armijos(@(y) ess1(x,y,lambda),a,pk);
    
    if alpha == 0
        lm = lm*beta;
        fprintf('backtracking failed to find the step, increasing the damping param: %d \n',lm);
        continue;
    else
        % update
        fprintf('step size=%g,',alpha);
        a = a+alpha*pk;
        lm = lm/beta;
    end
    stat(1,k) = r'*r;
    stat(2,k) = norm(J'*r);
    k=k+1;
end
figure; plot(stat(1,:)); title('objective values');
figure; plot(stat(2,:)); title('gradient norms');
end

function a = minsearch(x,K,a)
% User quadratic penalty method to solve the nonlinear constrained
% minimization
% f(x) = f(x) + .5*mu*(tx*tx+ty*ty+tz*tz-1)

x1 = h2e(K\[x(1:2,:);ones(1,size(x,2))]);
x2 = h2e(K\[x(3:4,:);ones(1,size(x,2))]);
x = [x1;x2];

mu = 1;
for j=1:20
    
    for k = 1:100
        [f,grad,H] = log_barrier(x,a,mu);
        if norm(grad) < 1/k
            break;
        end
        
        % search direction: steepest descent
        pk = -H\grad;
        
        % backtracking step size search
        alpha = armijos(@(y) log_barrier(x,y,mu),a,pk);
        
        % update
        a = a+alpha*pk;
    end
    stat(1,j) = f;
    stat(2,j) = norm(grad);
    mu = mu*10;
end
figure; semilogy(stat(1,:),'.b'); title('objective values');
figure; semilogy(stat(2,:),'*r'); title('gradient norms');
end

function [fval,grad] = quad_penalty(x,a,mu)

[fval,grad,H] = ess2(x,a);
fval = fval+.5*mu*(tx*tx+ty*ty+tz*tz-1)^2;
grad = grad+2*mu*(tx*tx+ty*ty+tz*tz-1)*[0 0 0 tx ty tz]';
end

function [fval,fgrad,fH] = log_barrier(x,a,mu)
% -log(tx*tx+ty*ty+tz*tz-1)
[fval,fgrad,fH] = ess2(x,a);

[tx,ty,tz] = deal(a(4),a(5),a(6));
cval = tx*tx+ty*ty+tz*tz-1;
cgrad = 2*[0 0 0 tx ty tz]';
cH = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];

bval = -log(cval);
bgrad = -cgrad/cval;
bH = cgrad*cgrad'/(cval*cval)-cH/cval;

fval = fval+mu*bval;
fgrad = fgrad+mu*bgrad;
fH = fH+bH;
end

function [fval,grad,H] = ess2(x,a)
% objective function for essential matrix estimation

% sum_j{r_j^2} where r_j = x_j'Ex_j and j runs over the correspondences

% grad r_j^2(k) = 2*r_j*drda_j
% hess r_j^2(i,k) = 2drda_i*drda_k+2r_j*dr2dra_idr_k

[rx,ry,rz,tx,ty,tz] = deal(a(1),a(2),a(3),a(4),a(5),a(6));
cx = cos(rx); sx = sin(rx);
cy = cos(ry); sy = sin(ry);
cz = cos(rz); sz = sin(rz);

w = [0 -tz ty; tz 0 -tx; -ty tx 0];
R = [cy*cz -cy*sz sy; sx*sy*cz+cx*sz -sx*sy*sz+cx*cz -sx*cy; -cx*sy*cz+sx*sz cx*sy*sz+sx*cz cx*cy];

Rdrx = [0 0 0; cx*sy*cz-sx*sz -cx*sy*sz-sx*cz -cx*cy; sx*sy*cz+cx*sz -sx*sy*sz+cx*cz -sx*cy];
Rdry = [-sy*cz sy*cz cy; cx*cy*cz -cx*cy*sz cx*sy; -sx*cy*cz sx*cy*sz -sx*sy];
Rdrz = [-cy*sz cy*sz 0; -sx*sy*sz+cx*cz -sx*sy*cz-cx*sz 0; cx*sy*sz+cx*sz cx*sy*cz-cx*sz 0];
wdtx = [0  0 0; 0 0 -1; 0 1 0];
wdty = [0  0 1; 0 0  0;-1 0 0];
wdtz = [0 -1 0; 1 0  0; 0 0 0];

Rdrxdrx = [0 0 0; -sx*sy*cz-cx*sz sx*sy*sz-cx*cz sx*cy; cx*sy*cz-sx*sz -sx*sy*sz-sx*cz -cx*sy];
Rdrxdry = [0 0 0; cx*cy*cz -cx*cy*sz cx*sy; sx*cy*cz -sx*cy*sz sx*sy];
Rdrxdrz = [0 0 0; -cx*sy*sz-sx*cz -cx*sy*cz+sx*sz 0; -sx*sy*sz+cx*cz -sx*sy*cz-cx*sz 0];

Rdrydrx = [0 0 0; -sx*cy*cz sx*cy*sz -sx*sy; -cx*cy*cz cx*cy*sz -cx*sy];
Rdrydry = [-cy*cz cy*cz -sy; -cx*sy*cz cx*sy*sz cx*cy; sx*sy*cz -sx*sy*sz sx*cy];
Rdrydrz = [sy*sz -sy*sz cy; -cx*cy*sz -cx*cy*cz 0; sx*cy*sz sx*cy*cz 0];

Rdrzdrx = [0 0 0; sx*sy*sz-cx*cz sx*sy*cz+cx*sz 0; -cx*sy*sz-sx*cz -cx*sy*cz+sx*sz 0];
Rdrzdry = [0 0 0; -cx*cy*sz -cx*sy*cz 0; -sx*cy*sz -sx*cy*cz 0];
Rdrzdrz = [0 0 0; -cx*sy*cz+sx*sz cx*sy*sz+sx*cz 0; -sx*sy*cz-cx*sz sx*sy*sz-cx*cz 0];

grad = zeros(6,1); H = zeros(6); fval = 0;
for j=1:size(x,2)
    y1 = [x(1:2,j);1]; y2 = [x(3:4,j);1];
    r = y1'*w*R*y2;
    
    % grad r_j
    rdrx = y1'*w*Rdrx*y2;
    rdry = y1'*w*Rdry*y2;
    rdrz = y1'*w*Rdrz*y2;
    rdtx = y1'*wdtx*R*y2;
    rdty = y1'*wdty*R*y2;
    rdtz = y1'*wdtz*R*y2;
    
    h11 = 2*rdrx*rdrx+2*r*y1'*w*Rdrxdrx*y2;
    h12 = 2*rdrx*rdry+2*r*y1'*w*Rdrxdry*y2;
    h13 = 2*rdrx*rdrz+2*r*y1'*w*Rdrxdrz*y2;
    h14 = 2*rdrx*rdtx+2*r*y1'*wdtx*Rdrx*y2;
    h15 = 2*rdrx*rdty+2*r*y1'*wdty*Rdrx*y2;
    h16 = 2*rdrx*rdtz+2*r*y1'*wdtz*Rdrx*y2;
    
    h21 = 2*rdry*rdrx+2*r*y1'*w*Rdrydrx*y2;
    h22 = 2*rdry*rdry+2*r*y1'*w*Rdrydry*y2;
    h23 = 2*rdry*rdrz+2*r*y1'*w*Rdrydrz*y2;
    h24 = 2*rdry*rdtx+2*r*y1'*wdtx*Rdry*y2;
    h25 = 2*rdry*rdty+2*r*y1'*wdty*Rdry*y2;
    h26 = 2*rdry*rdtz+2*r*y1'*wdtz*Rdry*y2;
    
    h31 = 2*rdrz*rdrx+2*r*y1'*w*Rdrzdrx*y2;
    h32 = 2*rdrz*rdry+2*r*y1'*w*Rdrzdry*y2;
    h33 = 2*rdrz*rdrz+2*r*y1'*w*Rdrzdrz*y2;
    h34 = 2*rdrz*rdtx+2*r*y1'*wdtx*Rdrz*y2;
    h35 = 2*rdrz*rdty+2*r*y1'*wdty*Rdrz*y2;
    h36 = 2*rdrz*rdtz+2*r*y1'*wdtz*Rdrz*y2;
    
    h41 = 2*rdtx*rdrx+2*r*y1'*wdtx*Rdrx*y2;
    h42 = 2*rdtx*rdry+2*r*y1'*wdtx*Rdry*y2;
    h43 = 2*rdtx*rdrz+2*r*y1'*wdtx*Rdrz*y2;
    h44 = 2*rdtx*rdtx;
    h45 = 2*rdtx*rdty;
    h46 = 2*rdtx*rdtz;
    
    h51 = 2*rdty*rdrx+2*r*y1'*wdty*Rdrx*y2;
    h52 = 2*rdty*rdry+2*r*y1'*wdty*Rdry*y2;
    h53 = 2*rdty*rdrz+2*r*y1'*wdty*Rdrz*y2;
    h54 = 2*rdty*rdtx;
    h55 = 2*rdty*rdty;
    h56 = 2*rdty*rdtz;
    
    h61 = 2*rdtz*rdrx+2*r*y1'*wdtz*Rdrx*y2;
    h62 = 2*rdtz*rdry+2*r*y1'*wdtz*Rdry*y2;
    h63 = 2*rdtz*rdrz+2*r*y1'*wdtz*Rdrz*y2;
    h64 = 2*rdtz*rdtx;
    h65 = 2*rdtz*rdty;
    h66 = 2*rdtz*rdtz;
    
    grad(1) = grad(1)+2*r*rdrx;
    grad(2) = grad(2)+2*r*rdry;
    grad(3) = grad(3)+2*r*rdrz;
    grad(4) = grad(4)+2*r*rdtx;
    grad(5) = grad(5)+2*r*rdty;
    grad(6) = grad(6)+2*r*rdtz;
    
    H = H + [h11 h12 h13 h14 h15 h16; h21 h22 h23 h24 h25 h26; h31 h32 h33 h34 h35 h36;
        h41 h42 h43 h44 h45 h46; h51 h52 h53 h54 h55 h56; h61 h62 h63 h64 h65 h66];
    fval = fval+r*r;
end

end

function [f,grad,H] = sym_tr_obj(x,a,lambda)

syms rx ry rz tx ty tz;

Rx = [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)];
Ry = [cos(ry) 0 sin(ry); 0 1 0; -sin(ry) 0 cos(ry)];
Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];

R = Rx*Ry*Rz;

w = [0 -tz ty; tz 0 -tx; -ty tx 0];

E = w*R;
f = 0;
for y=x
    y1 = [y(1:2);1];
    y2 = [y(3:4);1];
    f = f+(y1'*E*y2)^2;
end

f = f + (lambda*(tx*tx+ty*ty+tz*tz-1))^2;

if nargout>1
    grad = gradient(f, [rx, ry, rz, tx, ty, tz]);
    grad = subs(grad,rx,a(1));
    grad = subs(grad,ry,a(2));
    grad = subs(grad,rz,a(3));
    grad = subs(grad,tx,a(4));
    grad = subs(grad,ty,a(5));
    grad = subs(grad,tz,a(6));
    grad = double(grad);
    
    if nargout>2
        H = hessian(f,[rx, ry, rz, tx, ty, tz]);
        H = subs(H,rx,a(1));
        H = subs(H,ry,a(2));
        H = subs(H,rz,a(3));
        H = subs(H,tx,a(4));
        H = subs(H,ty,a(5));
        H = subs(H,tz,a(6));
        H = double(H);
    end
end

f = subs(f,rx,a(1));
f = subs(f,ry,a(2));
f = subs(f,rz,a(3));
f = subs(f,tx,a(4));
f = subs(f,ty,a(5));
f = subs(f,tz,a(6));
f = double(f);
end

function [r,J] = gaussnewton_tr_obj(x,a,lambda)
% objective function for essential matrix estimation

% sum_j{r_j^2} where r_j = x_j'Ex_j and j runs over the correspondences

% grad r_j^2(k) = 2*r_j*drda_j
% hess r_j^2(i,k) = 2drda_i*drda_k+2r_j*dr2dra_idr_k

validate = 0;

[rx,ry,rz,tx,ty,tz] = deal(a(1),a(2),a(3),a(4),a(5),a(6));
cx = cos(rx); sx = sin(rx);
cy = cos(ry); sy = sin(ry);
cz = cos(rz); sz = sin(rz);

if validate
    syms srx sry srz stx sty stz
    sRx = [1 0 0; 0 cos(srx) -sin(srx); 0 sin(srx) cos(srx)];
    sRy = [cos(sry) 0 sin(sry); 0 1 0; -sin(sry) 0 cos(sry)];
    sRz = [cos(srz) -sin(srz) 0; sin(srz) cos(srz) 0; 0 0 1];
    sR = sRx*sRy*sRz;
    sw = [0 -stz sty; stz 0 -stx; -sty stx 0];
end

w = [0 -tz ty; tz 0 -tx; -ty tx 0];
R = [cy*cz -cy*sz sy; sx*sy*cz+cx*sz -sx*sy*sz+cx*cz -sx*cy; -cx*sy*cz+sx*sz cx*sy*sz+sx*cz cx*cy];

if validate
    g = subs(sR,srx,rx);
    g = subs(g,sry,ry);
    g = subs(g,srz,rz);
    g = subs(g,stx,tx);
    g = subs(g,sty,ty);
    g = subs(g,stz,tz);
    assert(norm(double(g)-R)<1e-12);
end

Rdrx = [0 0 0; cx*sy*cz-sx*sz, -cx*sy*sz-sx*cz, -cx*cy; sx*sy*cz+cx*sz, -sx*sy*sz+cx*cz, -sx*cy];
Rdry = [-sy*cz, sy*sz, cy; sx*cy*cz, -sx*cy*sz, cx*sy; -cx*cy*cz, cx*cy*sz, -cx*sy];
Rdrz = [-cy*sz, -cy*cz, 0; -sx*sy*sz+cx*cz, -sx*sy*cz-cx*sz 0; cx*sy*sz+cx*sz, cx*sy*cz-sx*sz 0];
wdtx = [0  0 0; 0 0 -1; 0 1 0];
wdty = [0  0 1; 0 0  0;-1 0 0];
wdtz = [0 -1 0; 1 0  0; 0 0 0];

if validate
    E = sw*sR;
end

N = size(x,2);
if nargout>1
    J = nan(N+1,6);
end
r = nan(N+1,1);
for j=1:size(x,2)
    y1 = [x(1:2,j);1]; y2 = [x(3:4,j);1];
    r(j) = y1'*w*R*y2;
    rdrx = y1'*w*Rdrx*y2;
    rdry = y1'*w*Rdry*y2;
    rdrz = y1'*w*Rdrz*y2;
    rdtx = y1'*wdtx*R*y2;
    rdty = y1'*wdty*R*y2;
    rdtz = y1'*wdtz*R*y2;
    
    if validate
        g1 = gradient((y1'*E*y2)^2,[srx sry srz stx sty stz]);
        g = subs(g1,srx,rx);
        g = subs(g,sry,ry);
        g = subs(g,srz,rz);
        g = subs(g,stx,tx);
        g = subs(g,sty,ty);
        g = subs(g,stz,tz);
        g = double(g);
    end
    
    if nargout>1
        J(j,1) = 2*r(j)*rdrx;
        J(j,2) = 2*r(j)*rdry;
        J(j,3) = 2*r(j)*rdrz;
        J(j,4) = 2*r(j)*rdtx;
        J(j,5) = 2*r(j)*rdty;
        J(j,6) = 2*r(j)*rdtz;
    end
end

r(N+1) = lambda*(tx*tx+ty*ty+tz*tz-1);
if nargout>1
    J(N+1,1) = 0;
    J(N+1,2) = 0;
    J(N+1,3) = 0;
    J(N+1,4) = 4*r(N+1)*tx;
    J(N+1,5) = 4*r(N+1)*ty;
    J(N+1,6) = 4*r(N+1)*tz;
end
end

function [val,grad] = ess1(x,a,lambda)

[r,J] = gaussnewton_tr_obj(x,a,lambda);
val = r'*r;
grad = J'*r;
end

function dogleg(fn,a0,opt_param)
% trust region algorithm
% seeks to minimize m(p) = f + g'p + p'Bp
delta = 1;
k=1;
for j=1:N
    
    fprintf('dogleg iter: %d:',j);
    
    % evaluate the function
    [r,J] = gaussnewton_tr_obj(x,a,lambda);
    
    fprintf('obj: %g, gnorm: %g,',r'*r,norm(J'*r));
    
    if norm(J'*r) < 1e-8
        break;
    end
    
    
    % Gauss Newton Hessian approximation
    H = J'*J;
    
    % Steepest descent direction
    pk = J'*r;
    
    % Newton search direction
    pn = -H\(J'*r);
    
    
    if alpha == 0
        lm = lm*beta;
        fprintf('backtracking failed to find the step, increasing the damping param: %d \n',lm);
        continue;
    else
        % update
        fprintf('step size=%g,',alpha);
        a = a+alpha*pk;
        lm = lm/beta;
    end
    stat(1,k) = r'*r;
    stat(2,k) = norm(J'*r);
    k=k+1;
end
figure; plot(stat(1,:)); title('objective values');
figure; plot(stat(2,:)); title('gradient norms');
end

function [R,t] = EtoRt(E,x1,x2,K)
[rot,t] = EssentialMatrixToCameraMatrix(E);
[R,t,~] = SelectCorrectEssentialCameraMatrix(rot,t,x1,x2,K);
end

function check_params(R,t,R1,t1)

t = t/norm(t);
t1 = t1/norm(t1);

assert(norm(skew(t)*R-skew(t1)*R1,'fro')<1e-12);
end








