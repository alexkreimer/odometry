function R = karcher_mean(RR)
tol = 1e-6;
N = length(RR);
R = RR{1};

while 1
    v = nan(3,N);
    for i=1:N
        vr = vrrotmat2vec(R'*RR{i})';
        v(:,i) = vr(1:3)*vr(4);
    end
    v = mean(v,2);
    if norm(v)<tol
        return;
    end
    R = R*vrrotvec2mat([v/norm(v);norm(v)]);
end

end
