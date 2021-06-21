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