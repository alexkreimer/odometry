function test_matching()

dbstop if error;
close all;

KITTI_HOME = '/home/kreimer/KITTI/dataset';
KITTI_HOME = fullfile('F:', 'KITTI' , 'dataset');
DBG_DIR = fullfile('F:', 'debug');

image_dir  = fullfile(KITTI_HOME, 'sequences', '00');
poses_file = fullfile(KITTI_HOME, 'poses','00.txt');

[i1, i2] = read_kitti_images(image_dir, 1);
[i3, i4] = read_kitti_images(image_dir, 2);

det_param.corner_num = 2000;
det_param.quality = .0001;

c1 = detect(i1, det_param);
c3 = detect(i3, det_param);

ext_param.r = 3;
match_param.r = 200;

[feat1, valid1] = extract(i1, c1, ext_param);
[feat3, valid3] = extract(i3, c3, ext_param);

c1 = c1(:, valid1);
c3 = c3(:, valid3);

info = struct('c1', c1, 'c2', c3, 'f1', feat1, 'f2', feat3);
m = match(c1, c3, feat1, feat3, match_param);

plot_matches(i1, i3, c1, c3(:, m(1,:)), 'left 1', 'left 2');

plot_corners(i1, c1, 'left 1');
plot_corners(i3, c3, 'left 2');

end



end

