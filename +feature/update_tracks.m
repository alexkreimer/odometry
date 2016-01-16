function M = update_tracks(features, P)
% M is a cell-array;
% M{j} has all tracklets of length j, found in cur frame.

% A = M{j} is an j\times N array (N is the number of tracklets)
% Each column of A contains indices of features that comprise the
% tracklet. A(k,l) is the index of the feature in frame i-j+k


empty = 1;
m12 = features.m12;
m21 = features.m21;

valid = zeros(1, size(features.c1,2));
for i = 1:length(m12)
    ind = find(m21(1,:)==m12(2,i),1);
    if isempty(ind)
        continue;
    end
    if m12(1,i) == m21(2,ind)
        empty = empty+1;
        valid(m12(1,i)) = 1;
    end
end

if nargin<2
    % tracks of len 1
    M{1,1} = find(valid);    
    % index of tracks
    M{2,1} = 1:length(M{1,1});
    return;
else
end


% max track len in previous frame
max_len = size(P,2);

% current matches
cur_mt = features.mt;
found = zeros(1, size(features.c1,2));
for j = 1:max_len
    % tracks of len j in previous frame
    prev_tracks = P{1,j};
    prev_indices= P{2,j};

    cur_tracks  = nan(j+1, 0);
    cur_indices = nan(1, 0);
    empty = 1;
    for k = 1:size(prev_tracks,2);
        ind = find(cur_mt(2,:) == prev_tracks(j,k), 1);
        if isempty(ind)
            % feature was not found in cur frame
            continue;
        end
        cur_tracks(1:j, empty) = prev_tracks(:,k);
        cur_tracks(j+1, empty) = cur_mt(1, ind);
        if ~valid(cur_mt(1,ind))
            continue;
        end
        found(cur_mt(1, ind)) = 1;
        cur_indices(1, empty) = prev_indices(1,k);
        empty = empty + 1;
    end
    M{1,j+1} = cur_tracks;
    M{2,j+1} = cur_indices;
end
M{1,1} = find(xor(valid,found));

n = length(M{1,1});
m = max(P{2,1});
M{2,1} = m+1:(n+m);
end
