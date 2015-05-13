function assertT(T1,T2,varargin)

if size(T1,1) == 3
    T1 = [T1;0 0 0 1];
end

if size(T2,1) == 3
    T2 = [T2;0 0 0 1];
end

T1 = unit_t(T1);
T2 = unit_t(T2);

if nargin>2 && strcmp(varargin{1},'inv')
    assert(norm(eye(4)-T1*T2,'fro')<1e-6);
else
    assert(norm(eye(4)-T1\T2,'fro')<1e-6);
end
