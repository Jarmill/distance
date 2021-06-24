%(Height, Center-y, Radius) of circles

%Height: minimum y coordinate on circle (as -H)

%inner circle
% h = 0.2;
% h = 0.8;
% h = 1.8;
h = 2.5;

%outer circle
% H = 1;
H = 3;
dist = 0.5;

x_moon = moon_contour_base(h, H, dist);

figure(38)
clf
hold on

%moon contour
plot(x_moon(1, :), x_moon(2, :), 'm', 'LineWidth', 3)

pbaspect([diff(xlim), diff(ylim), 1])

