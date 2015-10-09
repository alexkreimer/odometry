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
