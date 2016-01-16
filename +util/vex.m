function v = vex(S)
%VEX Convert skew-symmetric matrix to vector
%
% V = VEX(S) is the vector (3x1) which has the skew-symmetric matrix S (3x3)
%
%           | 0   -vz  vy|
%           | vz   0  -vx|
%           |-vy   vx  0 |
%
% Notes::
% - This is the inverse of the function SKEW().
% - No checking is done to ensure that the matrix is actually skew-symmetric.
% - The function takes the mean of the two elements that correspond to each unique
%   element of the matrix, ie. vx = 0.5*(S(3,2)-S(2,3))
%
% See also SKEW.

v = 0.5*[S(3,2)-S(2,3); S(1,3)-S(3,1); S(2,1)-S(1,2)];
end
