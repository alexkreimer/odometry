classdef patch_extractor < handle
    properties
        center = [];
        descriptor = [];
        valid = [];
    end
    methods
        function obj = patch_extractor(im,center,patchr)
            if nargin>0
                % extract a square patch centered @ center of size
                % (2*patch_sz+1)x(2*patch_sz+1)
                % if the center is too close to the border, valid will be false and
                % patch=[]
                x1 = center(1,:)-patchr;
                x2 = center(1,:)+patchr;
                y1 = center(2,:)-patchr;
                y2 = center(2,:)+patchr;
                obj.valid = x1>1 & y1>1 & x2<size(im,2) & y2<size(im,1);
                obj.center = center(:,obj.valid);
                for j=[find(obj.valid);1:size(obj.center,2)]
                    x = x1(j(1)):x2(j(1));
                    y = y1(j(1)):y2(j(1));
                    patch = im(y,x);
                    obj.descriptor(:,j(2)) = double(patch(:));
                end
            end
        end
    end
end