function test_tracks(DATASET_DIR, sequence)

D = dir(fullfile(DATASET_DIR, 'sequences', sequence, 'image_0', '*.png'));
num = length(D(not([D.isdir])));

tracks = {};
grey2rgb = @(im) cat(3, im, im, im);
for i = 1:num
    [i1, i2] = read_images(fullfile(DATASET_DIR, 'sequences', sequence), i);
    i1 = grey2rgb(i1);
    i2 = grey2rgb(i2);
    tracks = feature.read_tracks(fullfile(DATASET_DIR, 'tracks', 'ku', sequence, sprintf('%06d.txt',i)), tracks);
    plot_tracks(i1,i2,i,tracks);
end
end

function plot_tracks(i1,i2,frame_num,tracks)

for i=1:size(tracks,2)
    x = tracks{1,i};
    plot_matches(i1,i2,frame_num,i,x);
end
end

function plot_matches(i1,i2,frame_num, track_len, x)

f1 = figure;
subplot(211); imshow(i1); hold on; colormap jet;
subplot(212); imshow(i2); hold on; colormap jet;

pt = [x(end-3:end-2,:); x(end-1:end,:)];
[~, ind] = sort(pt(1,:));

pt = pt(:, ind);

subplot(211);
scatter(pt(1,:), pt(2,:), 10, pt(1,:));
title(sprintf('left image %d, current matches',frame_num, track_len));

subplot(212);
scatter(pt(3,:), pt(4,:), 10, pt(1, :));
title(sprintf('right image %d, current matches',frame_num, track_len));

f2 = figure;
subplot(211); imshow(i1); hold on;
subplot(212); imshow(i2); hold on;

n = size(x,1);
subplot(211);
p1 = x(1:4:n,:);
p2 = x(2:4:n,:);
if size(p1,1) == 1
    plot(p1,p2,'.');
else
    plot(p1,p2,'-o');
end
title(sprintf('left image %d, tracks of len %d',frame_num, track_len));

subplot(212);
p1 = x(3:4:n,:);
p2 = x(4:4:n,:);
if size(p1,1) == 1
    plot(p1,p2,'.');
else
    plot(p1,p2,'-o');
end
title(sprintf('right image %d, tracks of len %d',frame_num, track_len));

waitforbuttonpress;
close(f1);
close(f2);
end
