function r = vtheta2r(vtheta)
angle = norm(vtheta);
if angle == 0
    axis = [0 1 0];
else
    axis  = vtheta(1:3)/angle;
end
r = [axis angle];
end