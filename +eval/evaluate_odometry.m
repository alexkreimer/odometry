function [rotation_errors, translation_errors] = evaluate_odometry(sequence,res,algs)

close all;

%set(groot,'defaultAxesColorOrder',[1 0 0;0 1 0;0 0 1],...
%    'defaultAxesLineStyleOrder','-|--|:')
  
DATA_ROOT  = '/home/kreimer/KITTI/';
have_ss  = 0;
gt.poses = util.read_poses(fullfile(DATA_ROOT, 'dataset', 'poses', [sequence, '.txt']));    
for j = 1:length(res)
    data = strrep(res{j},'results','data');
    path = fullfile(DATA_ROOT,data,sequence);
    d = dir(path);
    [x,y]=sort(datenum(char({d.date})));
    most_recent_file=char(d(y(end)).name);
    load(fullfile(path,most_recent_file));
    for i = 1:length(algs)
        alg = algs{i};
        field = [res{j},'_',algs{i}];
        disp(field);
        if strcmp(field(end-1:end),'ss')
            if have_ss
                fprintf('skip %s\n', field);
                continue;
            else
                have_ss = 1;
            end
        end
        info.(field).poses = util.read_poses(fullfile(DATA_ROOT, res{j}, alg, 'data', [sequence, '.txt']));
    end
end

fields = fieldnames(info);
num_alg= length(fields);

display = cell(1,length(fields));
for i = 1:length(fields)
    display{i} = strrep(fields{i},'_',' ');
    display{i} = strrep(display{i}, 'results', 'IO');
    if strcmp(display{i}(end-1:end),'ss')
        display{i} = 'SS';
    end
end

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

figure; hold on; title(sprintf('sequence %s: bird-view of the path', sequence));
set(gca, 'LineStyleOrder', {'-', ':', '--', '-.'}); % different line styles
% Call hold so that the first plot doesn't affect the properties you just set
hold all



for i=1:num_alg
    field = fields{i};
    poses = info.(field).poses(:,:,1:num_frames);
    x = poses(1,4,:);
    z = poses(3,4,:);
    plot(x(:),z(:),'DisplayName',display{i});
end
poses = gt.poses(:,:,1:num_frames);
x = poses(1,4,:);
z = poses(3,4,:);
plot(x(:),z(:),'DisplayName','Ground Truth');

clickableLegend(cat(2,display,'Ground Truth'));
%legend(gca,'show'); 
axis equal;
xlabel('x [m]');
ylabel('z [m]');
width = 16; % cm 
height = 12; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])
print('-dpng','-r500',sprintf('top_path_%s', sequence));

%figure;
%hold on; title('trajectory distances');
%for i=1:num_alg
%    field = fields{i};
%    plot(info.(field).dist(1:num_frames),'DisplayName',display{i});
%end
%clickableLegend(display);
%legend(gca,'show');

%lengths = [100,200,300,400,500,600,700,800];
lengths = [10 ,20 ,30 ,40 ,50, 60 ,70 ,80 ];
%lengths = [1,2,3];
step_size = 10;

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
            
            p1_gt = [gt.poses(:,:,first_frame); 0 0 0 1];
            p2_gt = [gt.poses(:,:,last_frame); 0 0 0 1];
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
rotation_errors = rotation_errors./repmat(lengths,[num_alg 1]);
translation_errors = translation_errors./len_count;
translation_errors = 100*translation_errors./repmat(lengths,[num_alg 1]);

figure; hold on; title(sprintf('sequence %s: rotation error', sequence));
set(gca, 'LineStyleOrder', {'-', ':', '--', '-.'}); % different line styles
% Call hold so that the first plot doesn't affect the properties you just set
hold all

for k = 1:num_alg
    plot(lengths,rotation_errors(k,:), 'DisplayName', [display{k}, ' ', num2str(mean(rotation_errors(k,~isnan(rotation_errors(k,:)))))]);
    leg{k} = [display{k}, ' ', num2str(mean(rotation_errors(k,~isnan(rotation_errors(k,:)))))];
end
ylabel('deg/m');
legend(gca,'show');

clickableLegend(leg);
width = 16; % cm 
height = 12; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])
print('-dpng','-r500',sprintf('rotation_error_%s', sequence));

figure; hold on; title(sprintf('sequence %s: translation error', sequence));
set(gca, 'LineStyleOrder', {'-', ':', '--', '-.'}); % different line styles
% Call hold so that the first plot doesn't affect the properties you just set
hold all
for k = 1:num_alg
    plot(lengths, translation_errors(k,:), 'DisplayName', [display{k}, ' ',num2str(mean(translation_errors(k,~isnan(translation_errors(k,:)))))]);
    leg{k} = [display{k}, ' ',num2str(mean(translation_errors(k,~isnan(translation_errors(k,:)))))];
end
ylabel('%');
clickableLegend(leg);
%legend(gca,'show');
width = 16; % cm 
height = 12; % cm
set(gcf, 'PaperPosition', [0, 0, width / 2.54, height / 2.54])
print('-dpng','-r500',sprintf('translation_error_%s', sequence));

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