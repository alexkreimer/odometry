function X = trg(x,param)

X = nan(3,size(x,2));

for i=1:size(x,2)
    X(:,i) = triangulate_naive(x(1:2,i),x(3:4,i),param.base,param.calib.f,param.calib.cu,param.calib.cv);
    if isnan(X(:,i))
        1;
    end
end

end
