function Rt = tr2mat(tr)
% XYZ

rx = tr(1); ry = tr(2); rz = tr(3);
tx = tr(4); ty = tr(5); tz = tr(6);

sx = sin(rx); cx = cos(rx); sy = sin(ry);
cy = cos(ry); sz = sin(rz); cz = cos(rz);

Rt(1,1) = +cy*cz;          Rt(1,2) = -cy*sz;          Rt(1,3) = +sy;    Rt(1,4) = tx;
Rt(2,1) = +sx*sy*cz+cx*sz; Rt(2,2) = -sx*sy*sz+cx*cz; Rt(2,3) = -sx*cy; Rt(2,4) = ty;
Rt(3,1) = -cx*sy*cz+sx*sz; Rt(3,2) = +cx*sy*sz+sx*cz; Rt(3,3) = +cx*cy; Rt(3,4) = tz;
Rt(4,1) = 0;               Rt(4,2) = 0;               Rt(4,3) = 0;      Rt(4,4) = 1;
