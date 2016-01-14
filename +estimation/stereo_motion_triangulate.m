function T = stereo_motion_triangulate(T1,T2,t0)
% combine motions
T1 = inv(T1);
T2 = inv(T2);

t1 = T1(1:3,4);
R1 = T1(1:3, 1:3);
t2 = T2(1:3,4);
R2 = T2(1:3,1:3);

P1 = [0 0 0]';
P2 = t0;
c = [-t1'*t1 t2'*t1; -t1'*t2 t2'*t2]\[ P2'*t1-P1'*t1; P2'*t2-P1'*t2];

P1 = P1 - c(1)*t1;
P2 = P2 - c(2)*t2;

t1 = P1 + (P2 - P1)/2;
t2 = t1 - t0;

% update motion direction of the cameras
T1 = [R1 t1; 0 0 0 1];
T = inv(T1);

end

function plot_epip(F, x1, x2, i1, i2, title1, title2)
figure;
subplot(211); imshow(i1); hold on;
title(title1);
plot(x1(1,:), x1(2,:),'go');
epiLines = epipolarLine(F, x2');
points = lineToBorderPoints(epiLines, size(i1));
line(points(:, [1,3])', points(:, [2,4])');
subplot(212); imshow(i2); hold on;
title(title2);
plot(x2(1,:), x2(2,:), 'go');
epiLines = epipolarLine(F', x1');
points = lineToBorderPoints(epiLines, size(i2));
line(points(:, [1,3])', points(:, [2,4])');
%truesize;
end


% symmetric signed distance from pts to the epilines
function d = dist_symm_epiline(F, x1, x2)

n1 = size(x1, 2);
n2 = size(x2, 2);
assert(n1 == n2);
d = nan(2, n1);
for i = 1:n1
    l2 = F*[x1(:, i); 1];
    l2 = l2/sqrt(l2(1)*l2(1)+l2(2)*l2(2));
    d1 = [x2(:, i); 1]'*l2;
    d(1, i) = d1;
    
    l1 = F'*[x2(:, i); 1];
    l1 = l1/sqrt(l1(1)*l1(1)+l1(2)*l1(2));
    
    d2 = [x1(:, i); 1]'*l1;
    d(2, i) = d2;
end

end



function txt = myupdatefcn(empt, event_obj, x, X)
% Customizes text of data tips

pos = get(event_obj,'Position');
d = x(3:4,:) - repmat(pos',[1 length(x)]);
[~, i] = min(sum(d.*d));
txt = {['Depth: ', num2str(X(3,i))], ['x: ', sprintf('%d %d', x(1,i), x(3,i))], ['d: ', sprintf('x:%d y:%d', x(1,i)-x(3,i), x(2,i)-x(4,i))]};
end
