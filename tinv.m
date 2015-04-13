function a1 = tinv(a)

a1 = nan(6,1);
T = inv(tr2mat(a));
[a1(1),a1(2),a1(3)] = decompose_rotation(T(1:3,1:3));
a1(4:6) = T(1:3,4);
end
