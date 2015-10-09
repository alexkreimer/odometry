function plot_corners(im, c, titl)
figure;
imshow(im);
hold on;
plot(c(1,:), c(2,:), '.g');
title(titl);
end

