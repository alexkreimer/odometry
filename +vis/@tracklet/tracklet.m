classdef tracklet < handle
    properties
        % array of feature objects, each object represents feature in a
        % different frame
        features;
        
        % length(obj.features)
        len = 0;
        
        % last time step when this tracklet was updated
        updated_at;
        
        % cache to speed things up?
        last_ft;
        
        % current time step
        step = 0;
        
        % matching tracklet in another view of a stereo pair
        match;
    end
    properties(Dependent)
        ppts;
        pts;
        descriptor;
        val;
    end
    methods
        function obj = tracklet(feature,step)
            if nargin>0
                if nargin<2
                    step = 1;
                end
                n = length(feature);
                obj(n) = vis.tracklet();
                for j=1:n
                    obj(j).features = feature(j);
                    obj(j).step = step;
                    obj(j).updated_at = obj.step;
                    obj(j).len = 1;
                    obj(j).last_ft = feature(j);
                end
            else
                obj.features = vis.feature.empty(1,0);
                obj.last_ft = vis.feature();
            end
        end
        function push_back(obj,feature)
            obj.features(end+1) = feature;
            obj.updated_at = obj.step;
            obj.len = length(obj.features);
            obj.last_ft = obj.features(obj.len);
        end
        
        function step_fwd(obj)
            n = numel(obj);
            if n == 1
                obj.step = obj.step+1;
            elseif n>1
                x = num2cell([obj.step]+1);
                [obj.step] = x{:};
            end
        end
        
        function feature = last(obj)
            % return the last seen feature of this tracklet
            if obj.len>0
                feature = obj.features(end);
            else
                feature = vis.feature();
            end
        end
        
        function feature = plast(obj)
            % returns one be4 last
            if obj.len>1
                feature = obj.features(obj.len-1);
            else
                feature = vis.feature();
            end
        end
        
        function v = valid(obj)
            % The tracklet is valid if it was updated at the current time
            % step, otherwise it was lost
            if numel(obj)>1
                v = [obj.updated_at] >= max([obj.step]);
            else
                v = obj.updated_at>=obj.step;
            end
        end
        
        function v = ppvalid(obj)
            % Checks if the tracklet is valid now (e.g., being tracked) and
            % has a previous value
            if numel(obj)>1
                v = and(obj.valid(),[obj.len]>1);
            else
                v = and(obj.valid(),obj.len>1);
            end
        end
        
        function update(obj,feature)
            % if the tracklet is valid and the feature is valid, update the
            % tracklet with the feature
            if and(obj.valid(t-1),feature.valid())
                obj.push_back(feature);
            end
        end
        
        function track_stupid(obj,im,searchr,patchr)
            % track features using simple template matching
            pts = [obj.pts];
            descriptors = [obj.descriptor];
            valid = [obj.valid];
            
            obj.step_fwd();
            parfor j=find(valid)
                % current feature point and its descriptor
                c = pts(:,j);
                d = descriptors(:,j);
                
                [x,y] = meshgrid(-searchr:+searchr,-searchr:+searchr);
                offset = [x(:),y(:)]';
                centers = bsxfun(@plus,offset,c);
                
                % extract all patches in a window
                extractor = vis.patch_extractor(im,centers,patchr);
                patches = extractor.descriptor;
                coords = extractor.center;
                
                % find min sad match
                [~,ind] = min(sum(abs(bsxfun(@minus,patches,d))));
                best_patch = patches(:,ind);
                best_pt = coords(:,ind);
                obj(j).push_back(vis.feature(best_pt, best_patch));
                
                %                 figure;
                %                 imshow(pim);
                %                 hold on;
                %                 plot(c(1),c(2),'or');
                %
                %                 figure;
                %                 imshow(im);
                %                 hold on;
                %                 plot(best_pt(1),best_pt(2),'og');
            end
        end
        
        function track(obj,im,searchr,pi1)
            pts = [obj.pts];
            valid = [obj.valid];
            
            descriptors = [obj.descriptor];
            obj.step_fwd();
            for j=find(valid)
                c = round(pts(:,j));
                d = descriptors(:,j);
                n = sqrt(length(d));
                template = reshape(d,[n n]);
                x1 = max(1,c(1)-searchr);
                x2 = min(size(im,2),c(1)+searchr);
                y1 = max(1,c(2)-searchr);
                y2 = min(size(im,1),c(2)+searchr);
                patch = im(y1:y2,x1:x2);
                [I_SSD,~] = template_matching(template,patch);
                [y,x] = find(I_SSD==max(I_SSD(:)));
                if length(x) > 1
                    % there is a number of matches, what should we do?
                    continue;
                end
                n = floor(n/2);
                
                %patch1 = patch(y-n:y+n,x-n:x+n);
                
                x = x1+x;
                y = y1+y;
                
                best_pt = [x;y];
                x1n = x-n-1;
                x2n = x+n-1;
                y1n = y-n-1;
                y2n = y+n-1;
                if x1n<1 || x2n>size(im,2) || y1n<1 || y2n>size(im,1)
                    continue;
                end
                best_patch = im(y1n:y2n,x1n:x2n);
                
                % debug
                %                 h = figure;
                %                 subplot(211); imshow(im); hold on; plot(x,y,'or'); rectangle('Position',[x1,y1,2*searchr+1,2*searchr+1]);
                %                 subplot(212); imshow(pi1); hold on; plot(c(1),c(2),'og');
                %                 figure;
                %                 subplot(241); imshow(patch);
                %                 subplot(242); imshow(template,[]);
                %                 subplot(243); imshow(best_patch,[]);
                %                 subplot(244); imshow(patch1,[]);
                %                 figure; imshow(I_SSD,[]);
                obj(j).push_back(vis.feature(best_pt, best_patch(:)));
            end
        end
        
        function pts = get.pts(obj)
            if obj.valid()
                pts = obj.last().pt;
            else
                pts = nan(2,1);
            end
        end
        
        function descriptor = get.descriptor(obj)
            if obj.valid()
                descriptor = obj.last().descriptor;
            else
                descriptor = nan(49,1);
            end
        end
        
        function val = get.val(obj)
            if obj.valid()
                val = true;
            else
                val = false;
            end
        end
        
        function pts = get.ppts(obj)
            if obj.ppvalid
                pts = [obj.last().pt;obj.plast().pt];
            else
                pts = nan(4,1);
            end
        end
        
        function augment(obj,im,num,margin,dist)
            % NOTE: this must be called for an array (e.g., obj is an array
            % of tracklets.  This method will search for new interest
            % points in image im and create a new tracklets for new
            % features.
            
            % detect features
            detector = vis.HarrisCorners(im,num,margin);
            
            % prune the points to remove the existing ones
            detector.prune([obj.pts],dist);
            
            % extract patches
            extractor = vis.patch_extractor(im,detector.corners,margin);
            
            % create new features
            new_features = vis.feature(extractor.center,extractor.descriptor);
            
            % update the tracklet array
            [n,m] = deal(length(obj),length(new_features));
            
            % preallocate
            obj(m+n) = vis.feature();
            
            % current time step
            cur_step = max([obj.step]);
            
            % append new tracklets
            for j=n+1:m
                obj(j) = vis.tracklet(feature,cur_step);
            end
        end
        
        function hmatch(obj,other,param)
            % match tracklets across the stereo pair
            
            % valid(live) tracklets for left cam
            valid1 = find([obj.valid]);
            % for the right cam
            valid2 = [other.valid];
            
            % collect image points
            pts1 = [obj.pts];
            pts2 = [other.pts];
            
            % for every (live) feature in the left image
            for i=1:length(valid1)
                min_ssd = inf;
                best_ft = nan;
                y = pts1(2,valid1(i));
                x = pts1(1,valid1(i));
                % find those that match epipolar constraint
                valid = find(abs(pts2(2,:)-y) <= param.threshy & abs(pts2(1,:)-x) <= param.threshx & valid2);
                ft1 = obj(valid1(i)).last;
                
                % find best (ssd) match
                for j=1:length(valid)
                    ft2 = other(valid(j)).last();
                    ssd = ft1.descriptor-ft2.descriptor;
                    ssd = sum(ssd.*ssd);
                    if ssd<min_ssd
                        best_ft = valid(j);
                        min_ssd = ssd;
                    end
                end
                
                if min_ssd<inf
                    % fill in the match information + 3D
                    f1 = obj(valid1(i)).last();
                    f2 = other(best_ft).last();
                    
                    if ~isempty(f2.match)
                        % f2 was already matched before
                        % it is possible to do better here
                        continue;
                    end
                    
                    f1.match = best_ft;
                    f2.match = valid1(i);

                    x1 = [x;y];
                    x2 = pts2(:,best_ft);
                    [X,d,~] = triangulate_naive(x1,x2,param.base,param.calib.f,param.calib.cu,param.calib.cv);
                    f1.X = X;
                    f1.disparity = d;
                    f2.X = X;
                    f2.disparity = d;
                end
            end
        end
        
        function [Xp,xp,x,ind] = get_matched(obj,other)
            % collect 3d data from previous frame and image points from
            % current frame to do estimate motion params
            k=1;
            for i=1:length(obj)
                if obj(i).valid() && length(obj(i).features)>1
                    % last feature in the track
                    f1 = obj(i).last();
                    
                    % one berfore last
                    pf1 = obj(i).plast();
                    
                    % check circle match & see if the disparity was not 0
                    if ~isempty(f1.match) && ~isempty(pf1.match) && ...
                            f1.match==pf1.match && other(f1.match).valid() && ...
                            ~any(isnan(pf1.X))
                        
                        % 3d in previous frame
                        Xp(:,k) = pf1.X;
                        
                        % image measuremens in current frame
                        x(:,k) = [f1.pt;other(f1.match).last().pt];
                        
                        % image meaurements in previous frame
                        xp(:,k) = [pf1.pt;other(f1.match).plast().pt];
                        
                        % index of this tracklet
                        ind(k) = i;
                        
                        k = k+1;
                    end
                end
            end
        end
        
        function pplot(obj,im,pim,titl)
            % show points that are valid now and are not new (e.g., were
            % tracked from the frame before)
            pts = [obj.ppts];
            valid = [obj.ppvalid];
            pts = pts(:,valid);
            showMatchedFeatures(im,pim,pts(1:2,:)',pts(3:4,:)','blend');
            title(titl);
        end
        
        function plot(obj,im,titl)
            % show points that are valid now and are not new (e.g., were
            % tracked from the frame before)
            pts = [obj.pts];
            valid = [obj.valid];
            pts = pts(:,valid);
            figure;
            imshow(im,[]);
            hold on;
            plot(pts(1,:),pts(2,:),'or');
            title(titl);
        end
    end
end
