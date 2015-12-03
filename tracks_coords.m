function [x,y] = tracks_coords(info, M, track_len, cur_frame)

if track_len>size(M,2)
    x = [];
    y = [];
    return;
end

all = [];
max_len = size(M,2);
for j = track_len:max_len
    tracks = M{j};
    first = j-track_len+1;
    last  = j;
    tracks = tracks(first:last,:);
    all = [all tracks];
end

n = size(all, 2);
x = nan(n, 1);
y = nan(n, 1);
for i = 1:n
    for j = 1:track_len
        f = all(j, i);
        c = info(cur_frame-track_len+j).c1(:,f);
        x(i) = c(1);
        y(i) = c(2);
    end
end
