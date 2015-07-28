function savePoses(filename, poses)
% save poses into txt file kitti style
fd = fopen(filename, 'w+');
for j=1:size(poses,3)
    fprintf(fd, '%g ', reshape(poses(1:3,:,j), [numel(poses(1:3,:,j)) 1]));
    fprintf(fd, '\n');
end
fclose(fd);
end