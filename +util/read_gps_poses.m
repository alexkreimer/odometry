function poses = read_gps_poses(file)

fd = fopen(file,'r');
poses = fscanf(fd,'%d %f %f %f', [4 Inf]);
fclose(fd);
poses = poses';
poses = poses(:,2:end);
p0 = repmat(poses(1,:),[size(poses,1) 1]);
poses = poses-p0;
end
