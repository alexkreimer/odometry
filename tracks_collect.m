function M = tracks_collect(features, cur_frame, P)

if nargin < 3
    n = length(features(cur_frame).c1);
    M{1} = 1:n;
    return;
end

max_len = size(P, 2);       % max tracklet len in previous frame
cur_mt = features(cur_frame).mt;    % current matches
found = zeros(1, size(features(cur_frame).c1,2));
for j = 1:max_len
    % tracklets of len j in previous frame
    prev_tracks = P{j};
    cur_tracks = nan(j+1, 0);
    empty = 1;
    for k = 1:size(prev_tracks, 2);
        ind = find(cur_mt(2,:) == prev_tracks(j,k), 1);
        if isempty(ind)
            % feature k was not found in cur frame
            continue;
        end
        cur_tracks(1:j, empty) = prev_tracks(:, k);
        cur_tracks(j+1, empty) = cur_mt(1, ind);
        found(cur_mt(1, ind)) = 1;
        empty = empty + 1;
    end
    M{j+1} = cur_tracks;
end
M{1} = find(~found);
end
