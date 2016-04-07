function harris_mono(ROOT_DIR, sequences, sha)

% * read image sequence
% * detect harris corners
% * match features

% DATASET_DIR = '/media/kreimer/my_drive/KITTI/dataset';

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
    
    D = dir([image_dir, '/image_0/*.png']);
    num = length(D(not([D.isdir])));

    for i = 1:num
        fprintf('sequence %s: processing frame %d of %d\n', sequence, i, num);
        i1 = read_images(image_dir, i-1);

        c1 = detect(i1, det_param);
        
        [feat1, valid1] = extract(i1, c1, ext_param);
        
        c1 = c1(:, valid1);

        info = struct('i1', i1, 'c1', c1, 'f1', feat1, 'm11p', []);

        x1 = info.c1;
        
        if ~first
            m11p = feature.match(c1, info_prv.c1, feat1, info_prv.f1, match_param_f);
            info.m11p = m11p;
        end
        first = false;
        tracks = feature.update_tracks_mono(info);
        feature.save_tracks(fullfile(out_dir,sprintf('%06d.txt',i-1)),tracks,x1,x2,m12,[]);
        info_prv = info;
        clear info;
    end
end
end