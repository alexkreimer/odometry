classdef FeatureDetector < handle
    % Interest point detector base class
    properties
        % 2xN numeric array that holds the image coordinates of the
        % detected points
        corners;
        
        % number of features to retain
        num;
        
        % features that are closer than margin to the boundary will be
        % discarded
        margin;
        
        % a copy of image
        im;
    end

    methods
        function obj = FeatureDetector(im,num,margin)
            if nargin>0
                obj.num = num;
                obj.margin = margin;
                obj.im = im;
                obj.corners = obj.detect();
            end
        end
        
        function valid = enforce_margin(obj,pts)
            % discard features that are closer than margin to the image
            % boundary
            row = pts(2,:);
            col = pts(1,:);
            [m,n] = size(obj.im);
            valid = row>obj.margin & row<m-obj.margin &...
                col>obj.margin & col<n-obj.margin;
        end

        function corners = filter_peaks(obj,filter)
            response = imfilter(obj.im,filter);
            bw = imregionalmax(response,8);
            bw = bwmorph(bw,'shrink',Inf);
            [row,col] = find(bw);
            corners = [col row]';
        end
    
        function prune(obj,pts,dist)
            % For every feature in obj.corners check if its dist pixels
            % close to any point in pts, if yes, it is considered invalid
            % and is removed
            n = size(obj.corners,2);
            valid = false(n,1);
            for j=1:n
                d = bsxfun(@minus,pts,obj.corners(:,j));
                if all(sum(d.*d)>dist*dist)
                    valid(j) = true;
                end
            end
            obj.corners = obj.corners(:,valid);
        end
    end
end