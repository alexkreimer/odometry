classdef horiz_match < handle
    properties
        match;
        X;
        d;
        epip_threshy = 2;
        epip_threshx = 50;
    end
    properties (Dependent)
        subs1;
        subs2;
        has3d;
    end
    methods
        function subs = get.subs1(obj)
            subs = obj.match(1,:);
        end
        function subs = get.subs2(obj)
            subs = obj.match(2,:);
        end
        function v = get.has3d(obj)
            v = isnan(obj.X(1,:));
        end
        function obj = horiz_match(i1,i2,ind,tracklets1,tracklets2,param)
            % when the constructor is done obj.match will be of size 2 x N
            % where N is the number of putative matches.  obj.match(1,:)
            % index into tracklets1, while obj.match(2,:) index into
            % tracklets2 (e.g., tracklets1(obj.match(1,i)) matches
            % tracklets2(obj.match(2,i))
            if nargin>0
                k=1;
                valid1 = find([tracklets1.valid]);
                valid2 = [tracklets2.valid];
                pts1 = [tracklets1.pts];
                pts2 = [tracklets2.pts];
                for i=1:length(valid1)
                    min_ssd = inf;
                    best_ft = nan;                    
                    y = pts1(2,valid1(i));
                    x = pts1(1,valid1(i));
                    valid = find(abs(pts2(2,:)-y)<=obj.epip_threshy &...
                        abs(pts2(1,:)-x)<=obj.epip_threshx &...
                        valid2);
                    ft1 = tracklets1(valid1(i)).last;
                    for j=1:length(valid)
                        ft2 = tracklets2(valid(j)).last();
                        ssd = ft1.descriptor-ft2.descriptor;
                        ssd = sum(ssd.*ssd);
                        if ssd<min_ssd
                            best_ft = valid(j);
                            min_ssd = ssd;
                        end
                    end
                    if min_ssd<inf
%                         h = figure;
%                         subplot(211); imshow(i1); hold on; plot(x,y,'og'); title('left'); hold off;
%                         subplot(212); imshow(i2); hold on; plot(pts2(1,best_ft),pts2(2,best_ft),'ob');hold off;
%                         saveas(h,sprintf('epip_match%02d.jpg',i)); close(h);
                        x1 = [x;y];
                        x2 = pts2(:,best_ft);
                        [X,d,~] = triangulate_naive(x1,x2,...
                            param.base,param.calib.f,param.calib.cu,param.calib.cv);
                        if d>param.min_disp
                            obj.X(:,k) = X;
                            obj.d(k) = d;
                            obj.match(1,k) = valid1(i);
                            obj.match(2,k) = best_ft;
                            
                            if isnan(tracklets1(valid1(i)).match)
                                % first time
                                tracklets1(valid1(i)) = best_ft;
                                tracklets2(best_ft) = valid1(i);
                            else
                                % check consistency with previous match
                                if tracklets1(valid1(i)) == best_ft
                                    % the match is consistent
                                end
                            end
                            
                            k=k+1;
                        end
                    end
                end
            end
        end
        
        function len = length(obj)
            len = size(obj.match,2);
        end
        
        function [ind,ind2] = get_ind2(obj,ind1)
            % if tracklets(ind1) <-> tracklets2(ind2) then given ind1 this
            % function returns ind2
            ind = find(obj.match(1,:)==ind1);
            if ~isempty(ind)
                ind2 = obj.match(2,ind);
            else
                ind2 = nan;
            end
        end
            
        function [cons,pcons] = get_cons(obj,prev)
            % find matches consistent with previous step
            len = obj.length();
            k = 1;
            for i=1:len
                ind = obj.match(1,i);
                [pind, pind2] = prev.get_ind2(ind);
                if ~isnan(pind2) && pind2 == obj.match(2,i)
                    cons(k) = i;
                    pcons(k) = pind;
                    k = k+1;
                end
            end
        end
        
        function new = copy(this)
            % Instantiate new object of the same class.
            new = feval(class(this));
            
            % Copy all non-hidden properties.
            p = properties(this);
            for i = 1:length(p)
                prop_info = findprop(this,p{i});
                if ~prop_info.Dependent
                    new.(p{i}) = this.(p{i});
                end
            end
        end
    end
end
