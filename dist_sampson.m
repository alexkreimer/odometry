function d = dist_sampson(F, x1, x2)
pts1h = e2h(x1);
pts2h = e2h(x2);

pfp = (pts2h' * F)';
pfp = pfp .* pts1h;
d = sum(pfp, 1) .^ 2;

epl1 = F * pts1h;
epl2 = F' * pts2h;
d = d ./ (epl1(1,:).^2 + epl1(2,:).^2 + epl2(1,:).^2 + epl2(2,:).^2);

end