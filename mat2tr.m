function a = mat2tr(T)
a = nan(6,1);
[a(1),a(2),a(3)] = decompose_rotation(T(1:3,1:3));
a(4:6) = T(1:3,4);
end