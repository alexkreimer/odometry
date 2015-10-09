function err = rot_error(pose_error)
a = pose_error(1,1);
b = pose_error(2,2);
c = pose_error(3,3);
d = 0.5*(a+b+c-1.0);
err = acos(max(min(d,1.0),-1.0));
end