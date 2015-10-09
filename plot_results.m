poses1 = kitti_read_poses('00_f.txt');
poses2 = kitti_read_poses('00_3d.txt');
poses3 = kitti_read_poses('F:\KITTI\dataset\poses\04.txt');

poses1 = permute(poses1(:,4,:), [1 3 2]);
poses2 = permute(poses2(:,4,:), [1 3 2]);

figure;
plot3(poses1(1,:), poses1(2,:), poses1(3,:), '.b');
hold on;
plot3(poses2(1,:), poses2(2,:), poses2(3,:), '.r');
hold on;
plot3(poses3(1,:), poses3(2,:), poses3(3,:), '.g');