function F = F_from_T(K, T)
t = T(1:3, 4);
R = T(1:3, 1:3);

F = inv(K')*skew(t)*R*inv(K);
end