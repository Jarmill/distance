% h_in = 1.9;
% h_out = 2;

% h_in = 2.6;
% h_in = 1;

% h_in = 2.99;
h_in_list = [0,0.3,0.5,0.8,1,1.2,1.5,2,2.5,2.99];

h_in = h_in_list(1);
h_out = 3;

%works when h_in, h_out <=1


figure(700)

for i = 1:length(h_in_list)
    h_in = h_in_list(i);
clf
hold on
%moon contour
x_moon = moon_base(h_in, h_out);

patch(x_moon(1, :), x_moon(2, :), 'k')
 
%% test the contour
c_in = [0.5*(1/h_in - h_in)];
r_in = 0.5*(1/h_in + h_in);

c_out = [0.5*(1/h_out - h_out)];
r_out = 0.5*(1/h_out + h_out);

dist = 1.05;

[x_inner, x_outer] = moon_intersection_points(c_in, c_out, dist);


x_cont = moon_contour_base_2(h_in, h_out, dist);

plot(x_cont(1, :), x_cont(2, :), 'r', 'LineWidth', 3)
pbaspect([diff(xlim), diff(ylim), 1])
axis off
title(sprintf('Inner Height = %0.2f', h_in), 'FontSize', 18)
pause

end