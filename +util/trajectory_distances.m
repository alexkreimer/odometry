function dist = trajectory_distances(poses)

N = size(poses,3);
dist = nan(N,1);
dist(1) = 0;
for i = 2:N
    p1 = poses(:,:,i-1);
    p2 = poses(:,:,i);
    
	dt = p2(1:3,4) - p1(1:3,4);
    
    dist(i) = dist(i-1)+sqrt(dt'*dt);
end

dist

end
