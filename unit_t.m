function T = unit_t(T)

T(1:3,4) = T(1:3,4)/norm(T(1:3,4));