function [match,info] = feat_match(i1,i2,f1,f2,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
addOptional(p,'do_sampson',false,@islogical);
addOptional(p,'sampson_thresh',2, @isnumeric);
addOptional(p,'F',eye(3));
addOptional(p,'norm_fn',@(x) sum(x.*x));
addOptional(p,'do_2nd_best', false, @islogical);
addOptional(p,'ratio', .6);
addOptional(p,'do_scale',false,@islogical);
addOptional(p,'scale_thresh',.05,@isnumeric);
addOptional(p,'do_orientation',false,@islogical);
addOptional(p,'orientation_thresh',cos(pi/360),@isnumeric);
addOptional(p,'debug',false,@islogical);

info.err_2nd_best = 0;
info.err_scale = 0;
info.err_orientation = 0;

parse(p,varargin{:});
match = nan(3,0);
for i = 1:length(f1)
    min_dst = inf;
    min_dst2= inf;
    min_ind = nan;
    if p.Results.do_sampson
        pt1 = f1(i).pt;
        x = [e2h(repmat(pt1, [1, length(f2)])); e2h([f2.pt])];
        [in,~] = funddist(p.Results.F, x, p.Results.sampson_thresh);
        if 0,
            figure;
            subplot(121); imshow(i1,[]); hold on; scatter(pt1(1),pt1(2),'red');
            subplot(122); imshow(i2,[]); hold on;
            for k=1:length(in)
                scatter(f2(in(k)).pt(1,:),f2(in(k)).pt(2,:));
            end
            epiLines = epipolarLine(p.Results.F, pt1');
            points = lineToBorderPoints(epiLines, size(i2));
            line(points(:, [1,3])', points(:, [2,4])');
            truesize;
        end
    else
        in = 1:length(f2);
    end
    if isempty(in),
        continue;
    end
    d1 = f1(i).d;
    for j = 1:length(in)
        d2 = f2(in(j)).d;
        dst = p.Results.norm_fn(d1-d2);
        if dst < min_dst,
            min_dst2 = min_dst;
            min_dst  = dst;
            min_ind  = in(j);
        end
    end
    if p.Results.do_2nd_best && min_dst > min_dst2*p.Results.ratio
        info.err_2nd_best = info.err_2nd_best + 1;
        continue;
    end
    if p.Results.do_scale
        v1 = reshape(f1(i).t,2,2)*[0;1];
        v2 = reshape(f2(in(j)).t,2,2)*[0;1];
        if abs(1.0-norm(v1)/norm(v2)) > p.Results.scale_thresh
            info.err_scale = info.err_scale + 1;
            continue;
        end
    end
    if p.Results.do_orientation
        v1 = reshape(f1(i).t,2,2)*[0;1];
        v2 = reshape(f2(in(j)).t,2,2)*[0;1];
        if abs(cos(0)-v1'*v2/(norm(v1)*norm(v2))) > p.Results.orientation_thresh
            info.err_orientation = info.err_orientation + 1;
            continue;
        end
    end

    match(1,end+1) = i;
    match(2,end) = min_ind;
    match(3,end) = min_dst;
end

% sort by the distance
[~,ind] = sort(match(3,:),'ascend');
match = match(:,ind);

if p.Results.debug
    figure;
    matchedPoints1 = round(p1(:,match(1,:)));
    matchedPoints2 = round(p2(:,match(2,:)));
    showMatchedFeatures(i1,i2,matchedPoints1',matchedPoints2','montage');
end
end

