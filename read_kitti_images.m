function [i1, i2] = read_kitti_images(seq_home, idx)
file_names = {fullfile(seq_home, 'image_0', sprintf('%06d.jpg',idx))
    fullfile(seq_home, 'image_1', sprintf('%06d.jpg',idx))};
i1 = imread(file_names{1});
i2 = imread(file_names{2});
end