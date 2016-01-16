function test_line_fit()

% epipole
e = [500; 180];

% number of lines
N = 5;

x = nan(2,N);
y = nan(2,N);
lines = nan(3,N);
lines_est = nan(3,N);
for i = 1:N
    % random line
    normal = randn(2,1);
    normal = normal/norm(normal);
    intersept = -normal(1)*e(1)-normal(2)*e(2);
    lines(:,i) = [normal; intersept];
    
    % generate data for current line
    x(:,i) = 500+500*randn(2,1);
    y(:,i) = get_line_points(normal, e, x(:,i));
end

a = [x(1,:); y(1,:)];
b = [x(2,:); y(2,:)];

figure;
hold on;
plot(e(1),e(2),'*r');
for i = 1:N
    x_min = min([a(1,i) b(1,i) e(1)]);
    x_max = max([a(1,i) b(1,i) e(1)]);
    xx = x_min:x_max;

    lines_est(:,i) = fit_line_e(a(:,i),b(:,i),e);
    yy = get_line_points(lines_est(1:2,i),e,xx);
    plot(xx,yy);
    plot(a(1,i), a(2,i),'og');
    plot(b(1,i), b(2,i),'ob');    
end
end
