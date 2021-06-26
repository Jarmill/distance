%(Height, Center-y, Radius) of circles

%Height: minimum y coordinate on circle (as -H)

%inner circle
% h = 0.9;
h = 0.8;
% h = 0.4;
% h = 1.3;
% h = 2.5;

%outer circle
% H = 3;
H = 2;
% H = 2;
% H = 3;

% dist = 6;
% dist = 0.5;

% [x_contour, x_moon] = moon_contour_base(h, H, dist);

figure(38)
clf

tl = tiledlayout(1,1);
ax = nexttile;
hold on

%moon contour
% plot(x_contour(1, :), x_contour(2, :), 'm', 'LineWidth', 3)
% plot(x_moon(1, :), x_moon(2, :), 'k', 'LineWidth', 3)



Ndist = 15;
dist_max = 4;

dist_range = linspace(0, dist_max, Ndist);
aut = autumn(length(dist_range));
x_moon = moon_base(h, H);

for i = 2:Ndist
%         x_contour = l1_circ_quart_contour(dist_range(i), R, 100);
    curr_dist = dist_range(i);

%          subplot(1,2, 1)   
    Xdist_curr = moon_contour_base(h, H, curr_dist);
    patch(Xdist_curr(1, :), Xdist_curr(2, :),'r', 'Linewidth', 3,'DisplayName', 'L1 Contour', ...        
    'EdgeColor', aut(i, :), 'FaceColor', 'None')
  

    
end


patch(x_moon(1, :), x_moon(2, :), 'r', 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')
pbaspect([diff(xlim), diff(ylim), 1])

set(ax, 'Colormap', autumn, 'CLim', [0, dist_max])
cbh = colorbar(ax(end)); 
% To position the colorbar as a global colorbar representing
% all tiles, 
cbh.Layout.Tile = 'east'; 
ylabel(cbh, 'Point-Set L1 Distance', 'FontSize', 16);

