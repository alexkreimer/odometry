function plot_triangles(i1, i1p, i2p, pt1, pt1p, pt2p)

pt1p(2,:) = pt1p(2,:) + size(i1,1);
pt2p(1,:) = pt2p(1,:) + size(i1,2);
pt2p(2,:) = pt2p(2,:) + size(i1,1);

rows = max(size(i1p,1), size(i2p,1));
cols = max(size(i1p,2), size(i2p,2));

im = zeros(2*rows, 2*cols);

[n,m] = size(i1);

im(1:n,1:m) = i1;
im(n+1:2*n, 1:m) = i1p;

[n1,m1] = size(i2p);
im(n+1:n+n1,m+1:m+m1) = i2p;

%im = [i1, 255*ones(size(i1)); i1p, i2p];
% imshow(im); hold on;
for i=1:size(pt1, 2)
    imshow(im,[]); hold on;
    plot([pt1(1,i) pt2p(1,i) pt1p(1,i) pt1(1,i)],...
        [pt1(2,i) pt2p(2,i) pt1p(2,i) pt1(2,i)]);
    hold off;
    waitforbuttonpress;
end
end
