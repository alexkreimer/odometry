function theta = angl(a,b)

a = a/norm(a);
b = b/norm(b);

theta = atan2(norm(cross(a,b)), dot(a,b));