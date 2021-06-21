%% contour plot of l1 distance of a circular arc
R = 1;          %radius of half-circle
% theta_c = pi/4;    %angular tilt of half-circle
% theta_c = 0;
% theta_c = pi;
% theta_c = 5*pi/6;
% theta_c = 5*pi/4;
% theta_c = 5*pi/3;
% theta_c = pi/3;
% theta_c = pi/6;

theta_list = [0; 5*pi/6; 5*pi/4; 5*pi/3];
% theta_c = 0.1;

Narc = 1000;
%unsafe set


% [Xdist, st_dist] = l1_halfcirc_contour_base(dist, R, theta_c, Narc);
% [Xdist, st_dist] = l1_halfcirc_contour(dist, R, theta_c, Narc);

%% plot the contour
figure(20)
clf

tl = tiledlayout(2, 2);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
% ax = cell(1, 4);
ax = [];

for j = 1:4
    ax(j) = nexttile;
hold on
theta_curr = theta_list(j);
%unsafe set
theta_half_range = linspace(theta_curr-pi/2, theta_curr + pi/2, Narc);
circ_half = [cos(theta_half_range); sin(theta_half_range)];
Xu = circ_half* R;
patch(Xu(1, :), Xu(2, :), 'r', 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')
% patch(Xdist(1, :), Xdist(2, :),'r', 'Linewidth', 3,'DisplayName', 'L1 Contour', 'EdgeColor', 'r', 'FaceColor', 'None')
    
%distance contours
Ndist = 15;
dist_max = 3;

dist_range = linspace(0, dist_max, Ndist);
aut = autumn(length(dist_range));

for i = 2:Ndist
%         x_contour = l1_circ_quart_contour(dist_range(i), R, 100);
    curr_dist = dist_range(i);

%          subplot(1,2, 1)   
    [Xdist_curr] = l1_halfcirc_contour(curr_dist, R, theta_curr, Narc);
    patch(Xdist_curr(1, :), Xdist_curr(2, :),'r', 'Linewidth', 3,'DisplayName', 'L1 Contour', ...        
    'EdgeColor', aut(i, :), 'FaceColor', 'None')
  
end

pbaspect([diff(xlim), diff(ylim), 1])
axis off
end
linkaxes(ax, 'xy');
set(ax, 'Colormap', autumn, 'CLim', [0, dist_max])
cbh = colorbar(ax(end)); 
% To position the colorbar as a global colorbar representing
% all tiles, 
cbh.Layout.Tile = 'east'; 
ylabel(cbh, 'Point-Set L1 Distance', 'FontSize', 16);
%% arc functions 

function [x_contour, st_contour]= l1_halfcirc_contour(dist, R, theta_c, Narc)
%     points with an l1 distance of dist to a half-circle with radius R
%     and angle theta_c (the midpoint of the circular arc is at theta_c)

    %handles transformations in case theta is outside the range [0, pi/4]

    %figure out the quadrant and remainder of the tilt angle
%     spins = floor(theta_c / (pi/2));
    
%     theta_c_mod = mod(theta_c, pi/2);
    
    theta_c_mod = wrapTo2Pi(theta_c);
    spins = 0;
    %different than modulus, due to the > sign in the while loop
    while theta_c_mod > pi/2
        theta_c_mod = theta_c_mod - (pi/2);
        spins = spins+1;
    end

    if theta_c_mod > pi/4
        theta_c_base = (pi/2) - theta_c_mod;
        [x_contour, st_contour] = l1_halfcirc_contour_base(dist, R, theta_c_base, Narc);
        x_contour = x_contour([2;1], :);
    else
        theta_c_base = theta_c_mod;
        [x_contour, st_contour] = l1_halfcirc_contour_base(dist, R, theta_c_base, Narc);        
    end
    
    %account for the rotation
    R_spin = [0, -1; 1, 0];
    x_contour = R_spin^(spins) * x_contour;
end

function [x_contour, st_contour] = l1_halfcirc_contour_base(dist, R, theta_c, Narc)
%     points with an l1 distance of dist to a half-circle with radius R
%     and angle theta_c (the midpoint of the circular arc is at theta_c)

    %assume that 0 <= theta_c <= pi/4

%     theta_c = pi/4;

    %right side arc
    Nright = Narc;
    theta_right = linspace(-pi/4, pi/4, Nright);
    arc_right = R*[cos(theta_right); sin(theta_right)] + [dist; 0];
    
    %top side arc
    Ntop = floor(Narc * (pi/4 + theta_c)/(pi/2))+1;
    theta_top = linspace(pi/4, pi/2 + theta_c, Ntop);
    arc_top = R*[cos(theta_top); sin(theta_top)] + [0; dist];
    
    %left side line
    line_left = [-R*sin(theta_c) - dist, R*sin(theta_c) - dist;
                  R*cos(theta_c), -R*cos(theta_c)];
              
    %bottom side arc
    Nbottom = floor(Narc * (pi/4 - theta_c) / (pi/2))+1;
    theta_bottom = linspace(-pi/4 -pi/4 + theta_c, -pi/4, Nbottom);
    arc_bottom = R*[cos(theta_bottom); sin(theta_bottom)] + [0; -dist];
    
    %bottom side arc meets right side arc (joined by line segment at corner
    %of half-circle)
    point_meet = [dist + sqrt(2)/2*R; -sqrt(2)/2*R];
    
    x_contour = [arc_right, arc_top, line_left, arc_bottom, point_meet]; 
    st_contour = struct('right', arc_right, 'top', arc_top, 'left', line_left, 'bottom', arc_bottom);
end