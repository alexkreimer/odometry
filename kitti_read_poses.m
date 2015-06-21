function poses = kitti_read_poses(poses_file)
fid = fopen(poses_file, 'r');
i=1;
tline = fgetl(fid);
while ischar(tline)
    poses(:,:,i) = reshape(sscanf(tline, '%f %f %f %f %f %f %f %f %f %f %f %f'),[4 3])';
    tline = fgetl(fid);
    i=i+1;
end
fclose(fid);
end
