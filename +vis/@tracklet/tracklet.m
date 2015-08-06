classdef tracklet < handle
    properties
        % array of feature objects, each object represents feature in a
        % different frame
        features;
        
        % length(obj.features)
        len = 0;
        
        % last time step when this tracklet was updated
        updated_at;
        created_at;
        
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
        function obj = tracklet(feature, step)
            if nargin>0
                if nargin<2
                    step = 1;
                end
                n = length(feature);
                obj(n) = vis.tracklet();
                for j=1:n
                    obj(j).features = feature(j);
                    obj(j).step = step;
                    obj(j).updated_at = step;
                    obj(j).created_at = step;
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

        function [features,n1] = tail(obj, n)
            % returns last n features if they exist
            % if there are less, n1 is the number of features returned
            if obj.len>n-1
                features = obj.features(1:n);
                n1 = n;
            else
                features = obj.features;
                n1 = obj.len;
            end
            return;
        end
        
        function v = valid(obj, step)
            % The tracklet is valid if it was updated at the current time
            % step, otherwise it was lost
            
            if nargin==1
                step = max([obj.step]);
            end
            
            v = [obj.updated_at] >= step & [obj.created_at] <= step;
        end
        
        function ft = get_feat(obj, step)
            ft = obj.features(step-obj.created_at+1);
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
        
        function track(obj, im, param, pi1)
            % track currently valid features in image 'im', search window
            % is given by searchr.  This will create a new feture in each
            % track that was successfully tracked

%             if param.do_dbg
%                 imshow(pi1);
%                 hold on;
%                 title('red are lost, green are tracked');
%             end
            pts = [obj.pts];
            valid = [obj.valid];
            
            fprintf('%d active tracks; ', numel(find(valid)));
            
            descriptors = [obj.descriptor];
            
            % this is important, since feature keep track of time
            obj.step_fwd();
            
            for j=find(valid)
                c = round(pts(:,j));
                
                % current descriptor
                d = descriptors(:,j);
                n = 2*param.patchr+1;
                template = reshape(d,[n n]);
                
                % search window in a new image
                x1 = max(1,c(1)-param.searchr);
                x2 = min(size(im,2),c(1)+param.searchr);
                y1 = max(1,c(2)-param.searchr);
                y2 = min(size(im,1),c(2)+param.searchr);
                patch = im(y1:y2,x1:x2);

                % do the search; normxcorr2 throws if the template is 'flat'
                try
                    cc = normxcorr2(double(template), double(patch));
                catch
                    if param.do_dbg
                        plot(c(1), c(2), '.r')
                    end
                    continue;
                end
                
                max_cc = max(cc(:));
                if max_cc < param.search_thresh
                    if param.do_dbg
                        plot(c(1), c(2), '.r')
                    end
                    
                    continue;
                else
                    [ypeak, xpeak] = find(cc==max_cc);
                end
                
                if length(xpeak)>1
                    if param.do_dbg
                        plot(c(1), c(2), '.r')
                    end
                    continue;
                end

                if param.do_dbg
                    plot(c(1), c(2), '.g')
                end

                xpeak = xpeak - param.patchr;
                ypeak = ypeak - param.patchr;

                n = param.patchr;
                x = x1 + xpeak;
                y = y1 + ypeak;

                best_pt = [x;y];
                x1n = x-n-1;
                x2n = x+n-1;
                y1n = y-n-1;
                y2n = y+n-1;
                if x1n<1 || x2n>size(im,2) || y1n<1 || y2n>size(im,1)
                    continue;
                end
                
                best_patch = im(y1n:y2n,x1n:x2n);
                obj(j).push_back(vis.feature(best_pt, best_patch(:)));
            end
            valid = [obj.valid];
            fprintf('%d successfully tracked\n', numel(find(valid)));
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
        
        function obj = augment(obj, im, num, margin, dist)
            % NOTE: this must be called for an array (e.g., obj is an array
            % of tracklets.  This method will search for new interest
            % points in image im and create a new tracklets for new
            % features.
            
            % detect features
            detector = vis.HarrisCorners(im, num, margin);
            
            % prune the points to remove the existing ones
            pts = obj.valid_pts();
            ptnum = size(pts,2);
            fprintf('%d valid tracks be4 augmentation; ', ptnum);
            detector.prune(pts, dist);
            
            % extract patches
            extractor = vis.patch_extractor(im,detector.corners, margin);
            
            % create new features
            new_features = vis.feature(extractor.center,extractor.descriptor);
            
            % update the tracklet array
            [n,m] = deal(length(obj),length(new_features));
            
            % preallocate
            if ptnum>2000 && numel(new_features)>10
                new_features = new_features(1:10);
                m = numel(new_features);
            end
            
            obj(m+n) = vis.tracklet();
            
            % current time step
            cur_step = max([obj.step]);
            
            % append new tracklets
            for j=1:m
                ft = vis.feature(new_features(j).frame, new_features(j).descriptor);
                obj(n+j) = vis.tracklet(ft, cur_step);
            end
            pts = obj.valid_pts();
            ptnum = size(pts,2);
            fprintf('%d after\n', ptnum);
        end
        
        function plot_matches(obj, other, i1, i2, titl)
            valid1 = find([obj.valid]);
            pts1 = [obj.pts];
            pts2 = [other.pts];
            k=1;
            for i = valid1
                m = obj(i).last().match;
                if ~isempty(m)
                    if isnan(other(m).match)
                        a = 1;
                    end
                    pt(1:2,k) = pts1(:,i);
                    pt(3:4,k) = pts2(:,obj(i).last().match);
                    k = k+1;
                end
            end
            c = linspace(1,255,length(pt));
            im1 = cat(3,i1,i1,i1);
            im2 = cat(3,i2,i2,i2);
            [~, ind] = sort(pt(1,:));
            pt(1,:) = pt(1,ind);
            pt(2,:) = pt(2,ind);
            pt(3,:) = pt(3,ind);
            pt(4,:) = pt(4,ind);
            figure;
            subplot(211); imshow(im1); hold on; scatter(pt(1, :), pt(2, :), 10, c); title(sprintf('%s; left', titl));
            subplot(212); imshow(im2); hold on; scatter(pt(3, :), pt(4, :), 10, c); title(sprintf('%s; right', titl));
        end
        
        function hmatch(obj, other, param)
            % match tracks across the stereo pair
            % after this function every valid feature has f.match set to a
            % matching feature in 'other' view
            %
            % TODO: try to check if 'symmetric' matching gives better
            % results
            
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
                        continue;
                    end
                    
                    f1.match = best_ft;
                    f2.match = valid1(i);

                    if isnan(f2.match)
                        a = 1;
                    end
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
        
        function [x, px, pX] = get_matched(obj, other, step)
            k=1;

            if nargin==2
                step = obj(1).step;
            end
            
            for i=1:length(obj)
                if obj(i).valid(step) && obj(i).valid(step-1)
                    % feature at time 'step' & 'step-1'
                    f1 = obj(i).get_feat(step);
                    pf1 = obj(i).get_feat(step-1);
                    
                    % check circle match & see if the disparity was not 0
                    if ~isempty(f1.match) && ~isempty(pf1.match) && ...
                            f1.match==pf1.match && other(f1.match).valid(step) && ...
                            ~any(isnan(pf1.X))
                        % 'other' features
                        other_f1 = other(f1.match).get_feat(step);
                        other_pf1= other(f1.match).get_feat(step-1);
                        
                        % 3d in previous frame
                        pX(:,k) = pf1.X;

                        % image measuremens in current frame
                        x(:,k) = [f1.pt; other_f1.pt];

                        % image meaurements in previous frame
                        px(:,k) = [pf1.pt; other_pf1.pt];

                        k = k+1;
                    end
                end
            end
        end
        
        function k = plot_circles(obj,i1,pi1,i2,pi2,other)
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

                        % image measuremens in current frame
                        x1(:,k) = f1.pt;
                        x2(:,k) = other(f1.match).last().pt;
                        px1(:,k) = pf1.pt;
                        px2(:,k) = other(f1.match).plast().pt;
                        k = k+1;
                    end
                end
            end
            
            px1(2,:) = px1(2,:) + size(i1,1);
            x2(1,:)  = x2(1,:)  + size(i1,2);
            px2(1,:) = px2(1,:) + size(i1,2);
            px2(2,:) = px2(2,:) + size(i1,1);
            
            im = [i1,i2;pi1,pi2];
            imshow(im); hold on;
            for i=1:k-1
                plot([x1(1,i) x2(1,i) px2(1,i) px1(1,i) x1(1,i)],...
                    [x1(2,i) x2(2,i) px2(2,i) px1(2,i) x1(2,i)]);
            end
            hold off;
        end
        function k = plot_circles1(obj,i1,pi1,i2,pi2,other)
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
                        
                        % image measuremens in current frame
                        x1(:,k) = f1.pt;
                        x2(:,k) = other(f1.match).last().pt;
                        px1(:,k) = pf1.pt;
                        px2(:,k) = other(f1.match).plast().pt;
                        k = k+1;
                    end
                end
            end
            
            [~, ind] = sort(x1(1,:));
            x1(1,:) = x1(1, ind);
            x1(2,:) = x1(2, ind);
            
            x2(1,:) = x2(1, ind);
            x2(2,:) = x2(2, ind);
            
            px1(1,:) = px1(1, ind);
            px1(2,:) = px1(2, ind);
            
            px2(1,:) = px2(1, ind);
            px2(2,:) = px2(2, ind);
            
            px1(2,:) = px1(2,:) + size(i1,1);
            x2(1,:)  = x2(1,:)  + size(i1,2);
            px2(1,:) = px2(1,:) + size(i1,2);
            px2(2,:) = px2(2,:) + size(i1,1);
            
            i1 = cat(3,i1,i1,i1);
            i2 = cat(3,i2,i2,i2);
            pi1 = cat(3,pi1,pi1,pi1);
            pi2 = cat(3,pi2,pi2,pi2);
            
            im = [i1,i2;pi1,pi2];
            c = linspace(1,255,length(x1));
            imshow(im); hold on;
            scatter(x1(1,:), x1(2,:), [], c);
            scatter(x2(1,:), x2(2,:), [], c);
            scatter(px1(1,:), px1(2,:), [], c);
            scatter(px2(1,:), px2(2,:), [], c);
            hold off;
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

        function pplot1(obj, im, prev_im, titl1, titl2)
            pts = [obj.ppts];
            valid = [obj.ppvalid];
            pts = pts(:,valid);
            
            im = cat(3,im,im,im);
            prev_im = cat(3,prev_im, prev_im, prev_im);
            
            [~, ind] = sort(pts(1,:));
            
            pts(1,:) = pts(1, ind);
            pts(2,:) = pts(2, ind);
            pts(3,:) = pts(3, ind);
            pts(4,:) = pts(4, ind);
            c = linspace(1, 255, length(pts));
            figure;
            subplot(211); imshow(im,[]); hold on; scatter(pts(1,:), pts(2,:), [], c); title(titl1);
            subplot(212); imshow(prev_im,[]); hold on; scatter(pts(3,:), pts(4,:), [], c); title(titl2)
            
            
        end

        function pts = valid_pts(obj)
            % get coords of valid features
            pts = [obj.pts];
            valid = [obj.valid];
            pts = pts(:,valid);
        end
        
        function plot(obj, im, titl)       
            % plot valid features
            pts = obj.valid_pts();
            imshow(im,[]);
            hold on;
            plot(pts(1,:), pts(2,:), '.g');
            title(titl);
        end
        
        function w = fit_line(obj,n)
            if nargin==1
                n = 3;
            end
            
            valid = [obj.valid];
            
            for i=find(valid)
                feat = obj(i).features;
                if length(feat)>=n
                    % last feature in the track
                    f1 = obj(i).last();
                    
                    % one berfore last
                    pf1 = obj(i).plast();
                    
                    % check circle match & see if the disparity was not 0
                    if ~isempty(f1.match) && ~isempty(pf1.match) && ...
                            f1.match==pf1.match
                        pts = [feat.pt];
                        % last 3 points
                        pts = pts(:,end-2:end);
                        x = pts(1,:)';
                        y = pts(2,:)';
                        A = [ones(3,1),x,y];
                        [c,n] = clsq(A,2);
                        plotline(x,y,'o',c,n,'--');
                        r = A*[c;n];
                        r_tot(i) = r'*r/numel(r);
                    end
                end
            end
            figure;
            plot(sort(r_tot));
        end
    end
end
