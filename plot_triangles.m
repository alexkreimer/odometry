function plot_triangles(i1, i1p, i2p, pt1, pt1p, pt2p)

pt1p(2,:) = pt1p(2,:) + size(i1,1);
pt2p(1,:) = pt2p(1,:) + size(i1,2);
pt2p(2,:) = pt2p(2,:) + size(i1,1);

im = [i1, 255*ones(size(i1)); i1p, i2p];
% imshow(im); hold on;
for i=1:size(pt1, 2)
    imshow(im); hold on;
    plot([pt1(1,i) pt2p(1,i) pt1p(1,i) pt1(1,i)],...
        [pt1(2,i) pt2p(2,i) pt1p(2,i) pt1(2,i)]);
    hold off;
    waitforbuttonpress;
end
end
