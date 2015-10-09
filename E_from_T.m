function E = E_from_T(T)
% T is a rigid transformation, E is the essential matrix

t = T(1:3, 4);
R = T(1:3, 1:3);
E = skew(t)*R;
E = E/E(3,3);

end