function [ratio, sigma] = cross_ratio(coords,e)

track_x = permute(coords(1,:,:),[2 3 1]);
track_y = permute(coords(2,:,:),[2 3 1]);

n  = size(track_x,2);
xr = nan(n,1);
yr = nan(n,1);
for i = 1:n
    [vx, vy] = deal(e(1),e(2));
    x  = track_x(1,i);
    xtt= track_x(2,i);
    xt = track_x(3,i);
    xr(i) = (xt-x)*(xtt-vx)/((xtt-xt)*(x-vx));
    
    y  = track_y(1,i);
    ytt= track_y(2,i);
    yt = track_y(3,i);
    yr(i) = (yt-y)*(ytt-vy)/((ytt-yt)*(y-vy));
end

if n>0
    todel = isinf(xr);
    xr(todel) = [];
    yr(todel) = [];
    todel = isinf(yr);
    xr(todel) = [];
    yr(todel) = [];
    
    % use both x corrdinate and y corrdinate ratios, same same
    data = [xr; yr];
    pd = fitdist(data, 'Normal');
    
    % plot the distribution
    %x_values = min(data):.1:max(data);
    %y = pdf(pd, x_values);
    %plot(x_values,y,'LineWidth',2)

    % save the result for output
    ratio = pd.mu;
    sigma = pd.sigma;
else
    ratio = nan;
    sigma = nan;
end
