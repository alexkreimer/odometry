function plot_matches(i1, i2, pt1, pt2, titl1, titl2)

grey2rgb = @(im) cat(3, im, im, im);
i1 = grey2rgb(i1);
i2 = grey2rgb(i2);

pt = [pt1; pt2];
[~, ind] = sort(pt(1,:));
pt = pt(:, ind);

figure;
subplot(211);
imshow(i1);
hold on;
scatter(pt(1,:), pt(2,:), 10, pt(1,:));
colormap jet;
title(titl1);

subplot(212);
imshow(i2);
hold on;
scatter(pt(3,:), pt(4,:), 10, pt(1, :));
colormap jet;
title(titl2);
end