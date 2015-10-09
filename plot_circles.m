function plot_circles(i1, i1_prv, i2, i2_prv, pt1, pt2, pt1_prv, pt2_prv)

pt1_prv(2,:) = pt1_prv(2,:) + size(i1,1);
pt2(1,:)  = pt2(1,:)  + size(i1,2);
pt2_prv(1,:) = pt2_prv(1,:) + size(i1,2);
pt2_prv(2,:) = pt2_prv(2,:) + size(i1,1);

im = [i1, i2; i1_prv, i2_prv];
% imshow(im); hold on;
for i=1:size(pt1, 2)
    imshow(im); hold on;
    plot([pt1(1,i) pt2(1,i) pt2_prv(1,i) pt1_prv(1,i) pt1(1,i)],...
        [pt1(2,i) pt2(2,i) pt2_prv(2,i) pt1_prv(2,i) pt1(2,i)]);
    hold off;
    waitforbuttonpress;
end
end
