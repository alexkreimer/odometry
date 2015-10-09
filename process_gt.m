function pgt = process_gt(poses_gt)
% KITTI gt holds pose of frame i as seen in frame 0

% Here we make them local, i.e. poses_gt(:, :, i) on output holds pose
% of frame {i-1} in frame {i}
num = size(poses_gt,3);
pgt = nan(3, 4, num);

% this would usually be I;0, i.e. the first frame location/orientation
% conincides with the global frame

pgt(:,:,1) = poses_gt(:,:,1);

for i=2:num
    T1 = [poses_gt(:,:,i); 0 0 0 1];      % frame {i} as seen from {0}
    T2 = [poses_gt(:,:,i-1); 0 0 0 1];    % frame {i-1} as seen from {0}
    T1 = inv(T1);                         % {0} as seen from {i}
    T = T1*T2;                            % frame {i-1} as seen from {i}
    pgt(1:3,:,i) = T(1:3,:);
end

end
