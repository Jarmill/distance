% h_in = 1.9;
% h_out = 2;

% h_in = 2.6;
% h_in = 1;

% h_in = 2.99;
% h_in_list = [0,0.3,0.5,0.8,1,1.2,1.5,2,2.5,2.99];
% 
% h_in = h_in_list(1);
% h_out = 3;

%circle heights of moon
h_in_list =  [0.4,1, 2.5]; 
h_out_list = [1, 2, 3];

%distance list
Ndist = 15;
dist_max = 3;

dist_range = linspace(0, dist_max, Ndist);
aut = autumn(length(dist_range));


%% plot the figure

figure(701)
clf
tl = tiledlayout(1, 3);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
ax = [];
for i = 1:3
    ax(i) = nexttile;

    h_in = h_in_list(i);
    h_out = h_out_list(i);
hold on
%moon contour
x_moon = moon_base(h_in, h_out);

patch(x_moon(1, :), x_moon(2, :), 'r', 'EdgeColor', 'None')
 
%% test the contour
% dist = 1.05;


for j = 2:Ndist
    dist_curr = dist_range(j);
x_cont = moon_contour_base_2(h_in, h_out, dist_curr);

% plot(x_cont(1, :), x_cont(2, :), 'r', 'LineWidth', 3)
%     [Xdist_curr] = l1_halfcirc_contour(curr_dist, R, theta_curr, Narc);
    patch(x_cont(1, :), x_cont(2, :),'r', 'Linewidth', 3,'DisplayName', 'L1 Contour', ...        
    'EdgeColor', aut(j, :), 'FaceColor', 'None')
end
pbaspect([diff(xlim), diff(ylim), 1])
axis off
% title(sprintf('Inner Height = %0.2f', h_in), 'FontSize', 18)
end


linkaxes(ax, 'xy');
set(ax, 'Colormap', autumn, 'CLim', [0, dist_max])
cbh = colorbar(ax(end)); 
% To position the colorbar as a global colorbar representing
% all tiles, 
cbh.Layout.Tile = 'east'; 
ylabel(cbh, 'Point-Set L2 Distance', 'FontSize', 16);