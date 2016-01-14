function y = get_line_points(n, v, x)
a = -n(1)/n(2);
b = n'*v/n(2);
y = a*x + b;
end