function [P0, P1] = kitti_read_calib(seq_home)
calib_file = fullfile(seq_home, 'calib.txt');
fd = fopen(calib_file, 'r');
p0 = fgetl(fd);
p1 = fgetl(fd);
P0 = reshape(sscanf(p0, 'P0: %f %f %f %f %f %f %f %f %f %f %f %f'),[4 3])';
P1 = reshape(sscanf(p1, 'P1: %f %f %f %f %f %f %f %f %f %f %f %f'),[4 3])';
fclose(fd);
end