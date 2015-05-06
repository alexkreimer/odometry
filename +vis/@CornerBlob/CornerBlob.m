classdef CornerBlob < vis.FeatureDetector
    % Interest point detector inspired by stereoscan.  It detects two kinds
    % of features, blobs and corners.
    properties (Access=private)
        % Geiger StereoScan masks
        blob_mask = [-1 -1 -1 -1 -1; -1 +1 +1 +1 -1; -1 +1 +8 +1 -1; -1 +1 +1 +1 -1; -1 -1 0 -1 -1];
        corn_mask = [-1 -1 +0 +1 +1; -1 -1  0 +1 +1;  0  0  0  0  0; +1 +1  0 -1 -1; +1 +1 0 -1 -1];
    end
    
    methods
        function obj = CornerBlob(im,num,margin)
            obj@vis.FeatureDetector(im,num,margin);
        end
        
        function corners = detect(obj)
            blobs = obj.filter_peaks(obj.blob_mask);
            corns = obj.filter_peaks(obj.corn_mask);
            
            if ~isempty(obj.margin)
                blobs = blobs(:,obj.enforce_margin(blobs));
                corns = corns(:,obj.enforce_margin(corns));
            end
            
            if ~isempty(obj.num)
                blobs = blobs(:,datasample(1:length(blobs),obj.num,'Replace',false));
                corns = corns(:,datasample(1:length(corns),obj.num,'Replace',false));
            end
            corners = [blobs corns];
        end
    end
end