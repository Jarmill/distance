h_in = 2;
h_out = 3;

%works when h_in, h_out <=1

x_moon = moon_base(h_in, h_out);

figure(700)
clf

%moon contour
plot(x_moon(1, :), x_moon(2, :), 'k', 'LineWidth', 3)
 

pbaspect([diff(xlim), diff(ylim), 1])