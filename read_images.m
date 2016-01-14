function [i1, i2] = read_images(seq_home, idx, file_mask)

if nargin < 3
    file_mask = '%06d.png';
end

cur_file = sprintf(file_mask,idx);

file1 = fullfile(seq_home, 'image_0', cur_file);
file2 = fullfile(seq_home, 'image_1', cur_file);

i1 = imread(file1);
i2 = imread(file2);
end