classdef tracklet < handle
    properties
        features;
        updated_at;
        last_ft;
        time = 0;
        len = 0;
    end
    methods
        function obj = tracklet
            obj.features = feat.empty(1,0);
            obj.last_ft = feat;
        end
        function push_back(obj,feature,t)
            obj.features(end+1) = feature;
            obj.updated_at = t;
            obj.len = obj.len + 1;
            obj.last_ft = obj.features(end);
        end
        function feature = last(obj)
            if obj.len>0
                feature = obj.features(end);
            else
                feature = feat;
            end
        end
        function feature = plast(obj)
            if obj.len>1
                feature = obj.features(obj.len-1);
            else
                feature = obj.last();
            end
        end
        function v = valid(obj,t)
            if nargin>1
                v = obj.updated_at>=t;
            else
                v = obj.updated_at>=obj.time;
            end
        end
        function len = length(obj)
            len = length([obj.features]);
        end
        function update(obj,feature,t)
            if t>obj.time
                obj.time = t;
            end
            if and(obj.valid(t-1),feature.valid())
                obj.push_back(feature,t);
            end
        end
    end
end
