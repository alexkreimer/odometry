classdef feat < handle
    properties
        frame;      %feature frame
        descriptor; %feature descriptor
        valid;
    end
    properties(Dependent)
        pt;
        x;
        y;
        d;
    end
    methods
        function obj = feat(frame,descriptor)
            if nargin>0
                obj.frame = reshape(frame,[2 1]);
                obj.descriptor = descriptor;
                obj.valid = true;
            else
                obj.valid = false;
            end
        end
        function pt = get.pt(obj)
            pt = nan(2,1);
            pt(1) = obj.x;
            pt(2) = obj.y;
        end
        function x = get.x(obj)
            if isempty(obj.frame)
                x = nan;
            else
                x = obj.frame(1);
            end
        end
        function y = get.y(obj)
            if isempty(obj.frame)
                y = nan;
            else
                y = obj.frame(2);
            end
        end
        function d = get.d(obj)
            d = obj.descriptor;
        end
    end
end
