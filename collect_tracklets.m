function [M, MX, MY] = collect_tracklets(info, num_frames)

num_tracks = size(info(2).mt,2);
M = zeros(num_frames,1000);

M (1, 1:num_tracks) = info(2).mt(2,:);
MX(1, 1:num_tracks) = info(1).c1(1, info(2).mt(2,:));
MY(1, 1:num_tracks) = info(1).c1(2, info(2).mt(2,:));

wm = num_tracks+1;
for i = 2:num_frames
    num_matches = size(info(i).mt,2);
    
    for j = 1:num_matches
        cur_ind = info(i).mt(1,j);
        prv_ind = info(i).mt(2,j);
        track_ind = find(M(i-1,:)==prv_ind);
        if ~isempty(track_ind)
            M(i, track_ind)  = cur_ind;
            MX(i, track_ind) = info(i).c1(1, cur_ind);
            MY(i, track_ind) = info(i).c1(2, cur_ind);
        else
            M(i,   wm) = cur_ind;
            M(i-1, wm) = prv_ind;
            MX(i,  wm) = info(i).c1(1, cur_ind);
            MY(i,  wm) = info(i).c1(2, cur_ind);
            MX(i-1,wm) = info(i-1).c1(1, prv_ind);
            MY(i-1,wm) = info(i-1).c1(2, prv_ind);
            wm = wm+1;
        end
    end
end
end
