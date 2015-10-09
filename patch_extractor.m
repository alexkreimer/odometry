function [patches, valid] = patch_extractor(im, c, ext_param)
% extract a square patch centered @ c of size (2*size+1)x(2*size+1)

r = ext_param.r;
x1 = c(1,:) - r;
x2 = c(1,:) + r;
y1 = c(2,:) - r;
y2 = c(2,:) + r;
valid = x1 > 1 & y1 > 1 & x2 < size(im,2) & y2 < size(im,1);
c = c(:, valid);
patches = nan((2*r+1)*(2*r+1), length(c));
for j = [find(valid); 1:length(c)]
    x = x1(j(1)):x2(j(1));
    y = y1(j(1)):y2(j(1));
    patch = im(y,x);
    patches(:, j(2)) = double(patch(:));
end
end

