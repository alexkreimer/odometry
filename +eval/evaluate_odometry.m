function evaluate_odometry(sequence)
close all;

DATA_ROOT  = '/media/kreimer/my_drive/KITTI/';

info.gt.poses = util.read_poses(fullfile(DATA_ROOT, 'dataset', 'poses', [sequence, '.txt']));
info.ss.poses = util.read_poses(fullfile(DATA_ROOT, 'results', 'ss', 'data', [sequence, '.txt']));
%info.no_opt_H.poses = util.read_poses(fullfile(DATA_ROOT, 'results', 'no_opt_H', 'data', [sequence, '.txt']));
%info.no_opt_F.poses = util.read_poses(fullfile(DATA_ROOT, 'results', 'no_opt_F', 'data', [sequence, '.txt']));
info.H_inf.poses = util.read_poses(fullfile(DATA_ROOT, 'results', 'H_inf', 'data', [sequence, '.txt']));

fields = fieldnames(info);
num_alg= length(fields);
num_frames = inf;
for i = 1:num_alg
    % there may be a different number of elements in each poses array
    field = fields{i};
    info.(field).dist = util.trajectory_distances(info.(field).poses);
    cur_frames = size(info.(field).poses,3);
    if num_frames > cur_frames
        num_frames = cur_frames;
    end
end

figure;
m=2;
n=2;

subplot(m,n,1); hold on; title('from-the top path view');
for i=1:num_alg
    field = fields{i};
    poses = info.(field).poses(:,:,1:num_frames);
    x = poses(1,4,:);
    z = poses(3,4,:);
    plot(x(:),z(:),'DisplayName',field);
end
legend(gca,'show'); axis equal;

subplot(m,n,2);
hold on; title('trajectory distances');
for i=1:num_alg
    field = fields{i};
    plot(info.(field).dist(1:num_frames),'DisplayName',field);
end
legend(gca,'show');

lengths = [100,200,300,400,500,600,700,800]';
lengths = [10,20,30,40,50 60,70,80];
lengths = [1,2,3];
step_size = 1;

rotation_errors = zeros(num_alg,length(lengths));
translation_errors = zeros(num_alg,length(lengths));

len_count = zeros(num_alg,length(lengths));

for l = 1:length(lengths)
    len = lengths(l);
    for first_frame = 1:step_size:num_frames
        for k = 1:num_alg
            field = fields{k};
            last_frame = get_last_frame_for_dist(info.(field).dist, first_frame, len);
            if isnan(last_frame)
                break;
            end
            len_count(k,l) = len_count(k,l) + 1;            
            
            p1_gt = [info.gt.poses(:,:,first_frame); 0 0 0 1];
            p2_gt = [info.gt.poses(:,:,last_frame); 0 0 0 1];
            pose_delta_gt = p1_gt\p2_gt;
            
            p1 = [info.(field).poses(:,:,first_frame); 0 0 0 1];
            p2 = [info.(field).poses(:,:,last_frame); 0 0 0 1];
            pose_delta = p1\p2;
            
            pose_error = pose_delta\pose_delta_gt;
            r_err = util.rot_error(pose_error);
            t_err = util.trans_error(pose_error);
            
            rotation_errors(k,l) = r_err + rotation_errors(k,l);
            translation_errors(k,l) = t_err + translation_errors(k,l);
        end
    end
end

rotation_errors = rotation_errors./len_count;
translation_errors = translation_errors./len_count;
translation_errors = translation_errors./repmat(lengths,[num_alg 1]);

subplot(m,n,3); hold on; title('rotation error');
for k = 1:num_alg
    field = fields{k};
    plot(lengths,rotation_errors(k,:), 'DisplayName', [strrep(field,'_',' '),' r error']);
end
legend(gca,'show');

subplot(m,n,4); hold on; title('translation error');
for k = 1:num_alg
    field = fields{k};
    plot(lengths, translation_errors(k,:), 'DisplayName', [strrep(field,'_',' '),' t error']);
end
legend(gca,'show');
end

function last_frame = get_last_frame_for_dist(dist, first_frame, len)

last_frame = nan;
for i=first_frame:length(dist)
    if dist(i)>dist(first_frame)+len
        last_frame = i;
        break;
    end
end
end