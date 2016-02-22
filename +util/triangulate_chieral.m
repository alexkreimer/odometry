function [X,visible] = triangulate_chieral(x1,x2,P1,P2,i1,i2)

num_pts = length(x1);
X = nan(4,num_pts);
visible = false(1,num_pts);
detM1 = det(P1(1:3,1:3));
detM2 = det(P2(1:3,1:3));
for j = 1:num_pts
    if nargin > 4
    figure;
    imshow([i1;i2]); hold on; plot([x1(1,j);x2(1,j)],[x1(2,j);x2(2,j)+size(i1,1)]);
    end
    
    %X(:, j) = vgg.vgg_X_from_xP_nonlin([x1(:, j) x2(:,j)],{P1,P2});
    X(:, j) = vgg.vgg_X_from_xP_lin([x1(:, j) x2(:,j)],{P1,P2});
    
    p1 = P1*X(:,j);
    p2 = P2*X(:,j);
    visible(j) = X(4,j)*p1(3)*detM1 > 0 && X(4,j)*p2(3)*detM2 > 0;
end
X = util.h2e(X);
end
