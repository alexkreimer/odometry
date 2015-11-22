function kitti_tracks_triangle()

% * read image sequence
% * detect harris corners
% * match features
% * prune outliers with triangle match heuristic

close all;
dbstop if error;

KITTI_HOME = '/home/kreimer/KITTI/dataset';
%KITTI_HOME = '/media/kreimer/my_drive/record_20150720/dataset';
DBG_DIR = fullfile('/home/kreimer/tmp/', 'debug');

sequence = '04';
image_dir  = fullfile(KITTI_HOME, 'sequences', sequence);

detect = @harris;
extract = @patch_extractor;

det_param.corner_num = 4000;
det_param.quality = .0001;
ext_param.r = 3;
match_param.r = 100;
hist = 100;

track_param.show_debug = 0;
track_param.save_debug = 0;

first = true;

% frame i-1:   1p ------ 2p
%              |         |
%              |         |
% frame i:     1 ------- 2

for i = 1:270
    fprintf('processing frame %d\n', i);
    [i1, i2] = read_kitti_images(image_dir, i);

    c1 = detect(i1, det_param);
    c2 = detect(i2, det_param);

    [feat1, valid1] = extract(i1, c1, ext_param);
    [feat2, valid2] = extract(i2, c2, ext_param);

    c1 = c1(:, valid1);
    c2 = c2(:, valid2);
    
    cur = i;
    prv = i-1;
    
    match_param.search1d = true;
    m12 = match(c1, c2, feat1, feat2, match_param);

    if track_param.show_debug
        for j = 1:size(m12, 2)
            showMatchedFeatures(i1, i2, c1(:, m12(1, j))', c2(:, m12(2, j))', 'montage');
            waitforbuttonpress;
        end
    end
    
    info(cur) = struct('i1', i1, 'i2', i2, 'c1', c1, 'c2', c2,...
        'f1', feat1, 'f2', feat2, 'm12', m12, 'm11p', [], 'm22p', [],...
        'm12p', [], 'm21p', [], 'mc', [], 'mt', []);
    
    if first == false
        match_param.search1d = false;
        
        % feature matching
        m11p = match(c1, info(prv).c1, feat1, info(prv).f1, match_param);
        m22p = match(c2, info(prv).c2, feat2, info(prv).f2, match_param);
        m21p = match(c2, info(prv).c1, feat2, info(prv).f1, match_param);
        m12p = match(c1, info(prv).c2, feat1, info(prv).f2, match_param);

        info(cur).m11p = m11p;
        info(cur).m22p = m22p;
        info(cur).m21p = m21p;
        info(cur).m12p = m12p;
        
        % pruning heuristics
        % circle
        mc = prune_circle(m12, m11p, m22p, info(prv).m12);
        
        % triangle
        mt = prune_triangle(info(prv).m12, m11p, m12p); % i1, i2, info(prv).i1), c1, c2, info(prv).c1);
        
        info(cur).mc = mc;
        info(cur).mt = mt;
        
%         pt1  = info(cur).c1(:, mt(1,:));
%         pt1p = info(prv).c1(:, mt(2,:));
%         pt2p = info(prv).c2(:, mt(3,:));
%         plot_triangles(i1, info(prv).i1, info(prv).i2, pt1, pt1p, pt2p);
        val = info(cur);
        save(['tracks/tracks_', sequence, '_', int2str(cur), '.mat'], 'val');
    else
        val = info(cur);
        save(['tracks/tracks_', sequence, '_', int2str(cur), '.mat'], 'val');        
        first = false;
    end
end
save(['tracks', sequence, '.mat'], 'info', '-v7.3');
end

function match = prune_circle(m12, m11, m22, m12_prv)

k = 1;
for i = 1:size(m12, 2)
    ind1 = m12(1, i);
    ind2 = m12(2, i);
    
    ind3 = find(m11(1,:) == ind1);
    ind4 = find(m22(1,:) == ind2);
    
    if or(isempty(ind3), isempty(ind4))
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

function match = prune_triangle(m12pp, m11p, m12p)
k = 1;
for i = 1:size(m12pp, 2)
    ind1 = m12pp(1, i);
    ind2 = m12pp(2, i);
    
    ind3 = find(m11p(2, :) == ind1);
    ind4 = find(m12p(2, :) == ind2);

    if or(isempty(ind3), isempty(ind4))
        continue
    end
    
    ind3 = m11p(1, ind3);
    ind4 = m12p(1, ind4);

    if ind3 == ind4
        match(:, k) = [ind3, ind1, ind2]';
        k = k + 1;
%         figure;
%         subplot(221); imshow(i1p); hold on; plot(c1p(1,ind3), c1p(2,ind3), 'or'); plot(c1p(1,ind4), c1p(2,ind4),'og');
%         subplot(223); imshow(i1); hold on; plot(c1(1,ind1), c1(2,ind1), 'og');
%         subplot(224); imshow(i2); hold on; plot(c2(1,ind2), c2(2,ind2), 'og');        
%         waitforbuttonpress;
    end
end

end