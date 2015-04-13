function [x,y,z] = decompose_rotation(R)
    % its assumed that R = R_z(z)R_y(y)R_x(x)
    
	%x = atan2(R(3,2), R(3,3));
	%y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
	%z = atan2(R(2,1), R(1,1));
    
    % http://en.wikipedia.org/wiki/Euler_angles#Intrinsic_rotations
    
    % R_x(x)R_y(y)R_z(z)
    x = atan2(-R(2,3), R(3,3));
    y = atan2(R(1,3), sqrt(R(1,1)*R(1,1)+R(1,2)*R(1,2)));
    z = atan2(-R(1,2),R(1,1));
end