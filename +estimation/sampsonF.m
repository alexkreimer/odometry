function val = sampsonF(F,x1,x2)

% x1, x2 s.t. x2'Fx1=0

N   = size(x1,2);
val = nan(N,1);
for i=1:N
    [p1, p2] = deal(x1(:,i), x2(:,i));

    p2tF = p2'*F;
    Fp1  = F*p1;
    v = sqrt(p2tF(1)^2 + p2tF(2)^2 + Fp1(1)^2  + Fp1(2)^2);
    if ~isnan(v)
        val(i) = p2'*F*p1/v;
    else
        val(i) = p2'*F*p1;
    end
end
end

