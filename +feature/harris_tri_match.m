function harris_tri_match(ROOT_DIR, sequences, sha)

% * read image sequence
% * detect harris corners
% * match features
% * prune outliers with triangle match heuristic

% e.g. DATASET_DIR = '/media/kreimer/my_drive/KITTI/dataset';

detect  = @feature.harris;
extract = @feature.patch_extractor;
    
det_param.corner_num = 4000;
det_param.quality = .0001;
ext_param.r = 3;
    
match_param_t.r = 100;
match_param_t.search1d = true;
    
match_param_f.r = 100;
match_param_f.search1d = false;
    
track_param.show_debug = 0;
    
for seq_num = 1:length(sequences)
    sequence  = sequences{seq_num};
    image_dir = fullfile(ROOT_DIR,'dataset','sequences', sequence);
    if nargin < 3
        out_dir = fullfile(ROOT_DIR, 'tracks', sequence);
    else
        out_dir = fullfile(ROOT_DIR, 'tracks', sha, sequence);
    end
    
    if exist(out_dir, 'dir')
        continue;
    else
        mkdir(out_dir);
    end
    
    first = true;
    
    % frame i-1:   1p ------ 2p
    %              |         |
    %              |         |
    % frame i:     1 ------- 2
    D = dir([image_dir, '/image_0/*.png']);
    num = length(D(not([D.isdir])));

    for i = 1:num
        fprintf('sequence %s: processing frame %d of %d\n', sequence, i, num);
        [i1, i2] = read_images(image_dir, i-1);
        
        c1 = detect(i1, det_param);
        c2 = detect(i2, det_param);
        
        [feat1, valid1] = extract(i1, c1, ext_param);
        [feat2, valid2] = extract(i2, c2, ext_param);
        
        c1 = c1(:, valid1);
        c2 = c2(:, valid2);
        m12 = feature.match(c1, c2, feat1, feat2, match_param_t);
        m21 = feature.match(c2, c1, feat2, feat1, match_param_t);

        if track_param.show_debug
            for j = 1:size(m12, 2)
                showMatchedFeatures(i1, i2, c1(:, m12(1, j))', c2(:, m12(2, j))', 'montage');
                waitforbuttonpress;
            end
        end
        
        info = struct('i1', i1, 'i2', i2, 'c1', c1, 'c2', c2,...
            'f1', feat1, 'f2', feat2, 'm12', m12, 'm21',m21,'m11p', [], 'm22p', [],...
            'm12p', [], 'm21p', [], 'mc', [], 'mt', []);

        x1 = info.c1;
        x2 = info.c2;
        
        if ~first
            % feature matching
            m11p = feature.match(c1, info_prv.c1, feat1, info_prv.f1, match_param_f);
            m22p = feature.match(c2, info_prv.c2, feat2, info_prv.f2, match_param_f);
            m21p = feature.match(c2, info_prv.c1, feat2, info_prv.f1, match_param_f);
            m12p = feature.match(c1, info_prv.c2, feat1, info_prv.f2, match_param_f);
            
            info.m11p = m11p;
            info.m22p = m22p;
            info.m21p = m21p;
            info.m12p = m12p;
            
            % pruning heuristics
            % circle
            %mc = prune_circle(m12, m11p, m22p, info_prv.m12);
            
            % triangle
            mt = prune_triangle(info_prv.m12, m11p, m12p, m12, m21);
            info.mt = mt;
            
%             pt1  = info.c1(:, mt(1,:));
%             pt2  = info.c2(:, mt(4,:));
%             pt1p = info_prv.c1(:, mt(2,:));
%             pt2p = info_prv.c2(:, mt(3,:));
%             util.plot_triangles(i1, i2, info_prv.i1, info_prv.i2, pt1, pt2, pt1p, pt2p);
            tracks = feature.update_tracks(info,tracks);
            
            feature.save_tracks(fullfile(out_dir,sprintf('%06d.txt',i-1)),tracks,...
                x1,x2,m12,mt);

            %save(fullfile(out_dir, ['frame_', int2str(i), '.mat']), 'info');
            info_prv = info;
            clear info;
        else
            %save(fullfile(out_dir, ['frame_', int2str(i), '.mat']), 'info');
            first = false;

            tracks = feature.update_tracks(info);
            feature.save_tracks(fullfile(out_dir,sprintf('%06d.txt',i-1)),tracks,...
                x1,x2,m12,[]);
            info_prv = info;
            clear info;
        end
    end
end

end

function match = prune_circle(m12, m11, m22, m12_prv)

match = nan(4,0);
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

function match = prune_triangle(m12pp, m11p, m12p, m12, m21)
% the order of values in mt : [cur_left; prv_left; prv_right]
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
        ind4 = find(m12(1,:) == ind3);
        ind5 = find(m21(2,:) == ind3);
        
        if isempty(ind4) || isempty(ind5)
            continue
        end
        
        if m12(2,ind4) ~= m21(1,ind5)
            continue;
        end
        match(:, k) = [ind3, ind1, ind2, m12(2,ind4)]';
        k = k + 1;
        %         figure;
        %         subplot(221); imshow(i1p); hold on; plot(c1p(1,ind3), c1p(2,ind3), 'or'); plot(c1p(1,ind4), c1p(2,ind4),'og');
        %         subplot(223); imshow(i1); hold on; plot(c1(1,ind1), c1(2,ind1), 'og');
        %         subplot(224); imshow(i2); hold on; plot(c2(1,ind2), c2(2,ind2), 'og');
        %         waitforbuttonpress;
    end
end

end