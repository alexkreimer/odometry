classdef feature < handle
    properties
        frame;      %feature frame(s)
        descriptor; %feature descriptor(s)
        X;
        disparity;
        match;      % index of a matching tracklet in a stereo pair
    end
    properties(Dependent)
        pt;
        x;
        y;
        d;
    end
    methods
        function obj = feature(frame,descriptor)
            switch nargin
                case 0
                    obj.frame = [];
                    obj.descriptor = [];
                case 1
                    n = size(frame,2);
                    if n==1
                        obj.frame = frame;
                    elseif n>1
                        obj(n) = vis.feature();
                        x = num2cell(frame,[1 n]);
                        [obj.frame] = x{:};
                    end
                case 2
                    n = size(frame,2);
                    m = size(descriptor,2);
                    assert(m==n);
                    if n==1
                        obj.frame = frame;
                        obj.descriptor = descriptor;
                    elseif n>1
                        obj(n) = vis.feature();
                        x = num2cell(frame,[1 n]);
                        [obj.frame] = x{:};
                        x = num2cell(descriptor,[1 n]);
                        [obj.descriptor] = x{:};
                    end
                otherwise
                    error('wrong number of args');
            end
        end
        
        function pt = get.pt(obj)
            if isempty(obj.frame)
                pt = [];
                return;
            end
            pt = obj.frame(1:2);
        end
        
        function x = get.x(obj)
            if isempty(obj.frame)
                x = [];
                return;
            end
            x = obj.frame(1);
        end
        function y = get.y(obj)
            if isempty(obj.frame)
                y = [];
                return;
            end
            y = obj.frame(2);
        end
        function d = get.d(obj)
            if isempty(obj.descriptor)
                d = [];
                return;
            end
            d = obj.descriptor;
        end
    end
end
