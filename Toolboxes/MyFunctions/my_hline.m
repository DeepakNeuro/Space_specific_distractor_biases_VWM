function r = my_hline(y_val, x_range)
% plot a horizontal line at y=y_val within x_range

g=ishold(gca);
hold on;

r = plot(x_range,[y_val, y_val]);
end