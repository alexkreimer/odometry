function test_est_F(pin)

x1 = pin.x1;
x2 = pin.x2;
inliers = pin.inliers;
n = length(pin.F);

h1 = 1;
h2 = 2;

h2 = figure(h2); 
subplot(211);
y = randsample(length(inliers), 10);
showMatchedFeatures(pin.i1, pin.i2, x1(:, inliers(y))', x2(:,inliers(y))', 'montage');
title(sprintf('Matches: frame (%d) %s vs frame (%d) %s', pin.ind1, pin.i1_name, pin.ind2, pin.i2_name));
legend({'1st', '2nd', 'match'});
subplot(212);

imshow(pin.i1);

if pin.dbg_save
    save_dbg(fullfile(pin.DBG_DIR, sprintf('matches_%04d_%04d.png', pin.ind1, pin.ind2)));
end

for j = 1:n
    F = pin.F{j};
    E = pin.E{j};
    
    fprintf('analysis of %s:\n', pin.name{j}');
    lambda = eig(E);
    lambda = lambda/lambda(2);
    fprintf('eigenvalues of E: %f, %f, %f\n', lambda(1), lambda(2), lambda(3));
    
    e = null(F);
    e = e/e(3);
    T = pin.T{j};
    R = T(1:3,1:3);
    t = T(1:3,4);
    e1 = pin.K*R'*t;
    e = e1/e1(3);
    fprintf('epipole of F: (%f, %f)\n', e(1), e(2));
    figure(h2); subplot(212); hold on; plot(e(1), e(2),'x'); hold off;
    
    [err, d] = residual_error(F, x1, x2, inliers);
    fprintf('residual error %g [px] for %d inliers out of %d points\n',...
        err, length(pin.inliers), length(pin.x1));
    
    h1 = figure(h1);
    subplot(1, n, j);
    width1 = .5;
    [N, edges] = histcounts(d(1, :));
    bar(edges(1:end-1), N, width1, 'FaceColor', [0.2,0.2,0.5]);
    
    if size(d,1) > 1
        [N, edges] = histcounts(d(2,:));
        bar(edges(1:end-1), N, width1/2, 'FaceColor',[0,0.7,0.7],'EdgeColor',[0,0.7,0.7]);
        title(sprintf('signed distance to the epipolar lines\n average %f', err));
    else
        title({sprintf('%s sampson distance', pin.name{j}), sprintf('average %f', err)});
    end

%     plot_epip(F, x1(:, inliers), x2(:, inliers), pin.i1, pin.i2, pin.name, '');
%     if pin.dbg_save
%         save_dbg(fullfile(pin.DBG_DIR, sprintf('epip_%04d_est_%d.png', pin.ind1,...
%             pin.ind2)));
%         close;
%     end
end

figure(h2);
subplot(212);
legend(pin.name);
end

function [err, d] = residual_error(F, x1, x2, inliers)
% computes average error from putative matches to the correspnding epipolar
% lines

dist_fn = @dist_sampson;

d = dist_fn(F, x1(:, inliers), x2(:, inliers));

err = sum(abs(d(:)))/length(inliers);

end

