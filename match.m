function m = match(c1, c2, f1, f2, match_param)

k = 1;
m = nan(3, 0);
for i = 1:length(c1)
    pt = c1(:, i);
    
    d = c2 - repmat(pt, [1 length(c2)]);
    d2 = sum(d.*d);
    
    if match_param.search1d
        % match of rectified images
        valid = find(d2 < match_param.r*match_param.r & abs(d(2,:)) < 2);
    else
        valid = find(d2 < match_param.r*match_param.r);
    end
    
    if isempty(valid)
        continue
    end
    
    ft = f1(:, i);
    d = f2(:, valid) - repmat(ft, [1 length(valid)]);
    d = sum(d.*d);
    [val, ind] = min(d);
    
    ind1 = find(m(2, :) == valid(ind));
    if ~isempty(ind1)
        if m(3, ind1) > val
            m(:, ind1) = [];
            k = k-1;
        else
            continue;
        end
    end
    m(1, k) = i;
    m(2, k) = valid(ind);
    m(3, k) = val;
    k = k + 1;
end