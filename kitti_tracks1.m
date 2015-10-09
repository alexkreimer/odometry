function kitti_tracks1()

% * read image sequence
% * detect harris corners
% * match features
% * prune outliers with circular match heuristic

close all;
dbstop if error;

KITTI_HOME = '/home/kreimer/KITTI/dataset';
KITTI_HOME = fullfile('F:', 'KITTI' , 'dataset');
DBG_DIR = fullfile('F:', 'debug');

image_dir  = fullfile(KITTI_HOME, 'sequences', '00');
poses_file = fullfile(KITTI_HOME, 'poses','00.txt');

detect = @harris;
extract = @patch_extractor;

det_param.corner_num = 2000;
det_param.quality = .0001;
ext_param.r = 3;
match_param.r = 100;
hist = 2;

track_param.show_debug = 1;
track_param.save_debug = 0;

first = true;

for i = 1:250
    fprintf('processing frame %d\n', i);
    [i1, i2] = read_kitti_images(image_dir, i);

    c1 = detect(i1, det_param);
    c2 = detect(i2, det_param);

    [feat1, valid1] = extract(i1, c1, ext_param);
    [feat2, valid2] = extract(i2, c2, ext_param);

    c1 = c1(:, valid1);
    c2 = c2(:, valid2);
    
    cur = mod(i-1, hist) + 1;
    prv = mod(i-2, hist) + 1;
    
    match_param.search1d = true;
    m12 = match(c1, c2, feat1, feat2, match_param);

    if track_param.show_debug
        for j = 1:size(m12, 2)
            showMatchedFeatures(i1, i2, c1(:, m12(1, j))', c2(:, m12(2, j))', 'montage');
            waitforbuttonpress;
        end
    end
    
    info(cur) = struct('i1', i1, 'i2', i2, 'c1', c1, 'c2', c2,...
        'f1', feat1, 'f2', feat2, 'm12', m12, 'm11', [], 'm22', []);
    
    if first == false
        match_param.search1d = false;
        m11 = match(c1, info(prv).c1, feat1, info(prv).f1, match_param);
        
        if param.show_debug
            figure; 
            plot(sort(m11(3,:)));
            title('left vs prev left: L2 distances betwen matched patches')
        end

         for j = 1:size(m11, 2)
             showMatchedFeatures(i1, info(prv).i1, c1(:, m11(1, j))', info(prv).c1(:, m11(2, j))');
             title('left vs. left previous')
             waitforbuttonpress;
         end
        
        m22 = match(c2, info(prv).c2, feat2, info(prv).f2, match_param);
        
        if param.show_debug
            figure;
            plot(sort(m22(3,:)));
            title('right vs prev right: L2 distances between matched patches');
        end
        
        info(cur).m11 = m11;
        info(cur).m22 = m22;
        
        mc = match_circle(m12, m11, m22, info(prv).m12);
        
        pt1 = info(cur).c1(:, mc(1,:));
        pt2 = info(cur).c2(:, mc(2,:));
        pt1_prv = info(prv).c1(:, mc(3,:));
        pt2_prv = info(prv).c2(:, mc(4,:));
%         plot_circles(i1, info(prv).i1, i2, info(prv).i2, pt1, pt2, pt1_prv, pt2_prv);
    else
        first = false;
    end
end

end

function match = match_circle(m12, m11, m22, m12_prv)

k = 1;
for i = 1:size(m12, 2)
    ind1 = m12(1, i);
    ind2 = m12(2, i);
    
    ind3 = find(m11(1,:) == ind1);
    ind4 = find(m22(1,:) == ind2);
    
    if or(isempty(ind3), isempty(ind2))
        continue
    end
    
    ind3 = m11(2, ind3);
    ind4 = m22(2, ind4);
    
    ind = find(m12_prv(1,:) == ind3);
    
    if ~isempty(ind) && m12_prv(2, ind) == ind4
        match(:, k) = [ind1, ind2, ind3, ind4]';
        k = k + 1;
    end
end

end