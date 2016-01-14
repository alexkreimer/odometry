function savePoses(filename, poses)
% save poses into txt file kitti style
fd = fopen(filename, 'w+');
for i=1:size(poses,3)
    fprintf(fd, '%g %g %g %g %g %g %g %g %g %g %g %g\n', poses(1,1,i), poses(1,2,i), poses(1,3,i), poses(1,4,i),...
    poses(2,1,i), poses(2,2,i), poses(2,3,i), poses(2,4,i),...
    poses(3,1,i), poses(3,2,i), poses(3,3,i), poses(3,4,i));
end
fclose(fd);
end