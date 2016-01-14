function poses = read_poses(poses_file)
fid = fopen(poses_file, 'r');
i=1;
tline = fgetl(fid);
while ischar(tline)
    val = sscanf(tline, '%f %f %f %f %f %f %f %f %f %f %f %f');
    poses(1,1,i) = val(1);
    poses(1,2,i) = val(2);
    poses(1,3,i) = val(3);
    poses(1,4,i) = val(4);

    poses(2,1,i) = val(5);
    poses(2,2,i) = val(6);
    poses(2,3,i) = val(7);
    poses(2,4,i) = val(8);

    poses(3,1,i) = val(9);
    poses(3,2,i) = val(10);
    poses(3,3,i) = val(11);
    poses(3,4,i) = val(12);
    
    tline = fgetl(fid);
    i=i+1;
end
fclose(fid);
end
