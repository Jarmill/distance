%quarter-torus
R = 1;  %radius of half-sphere
r = 1; %inner radius of torus

Nth = 30;
Nphi = 40;

% phi = linspace(4*pi/2, 3*pi, Nphi);
phi = linspace(0, 2*pi, Nphi);  %perimiter
% theta = linspace(0, 2*pi, Nth); %azimuth
theta = linspace(0, pi/2, Nth);

[Phi, Theta] = meshgrid(phi, theta);

%
X_torus = (R+r.*cos(Theta)).*cos(Phi);
Y_torus = (R+r.*cos(Theta)).*sin(Phi);
Z_torus = r.*sin(Theta);

X_sphere = (R+r)*cos(Theta).*cos(Phi);
Y_sphere = (R+r)*cos(Theta).*sin(Phi);
Z_sphere = -(R+r)*sin(Theta);

rscale = R/(R+r);
circ = R*[cos(phi); sin(phi)];
z_circ = r * ones(size(phi));
% [X_torus, Y_torus, Z_torus] = pol2cart(Theta, Phi, Z_top);

figure(1)
clf

hold on
falpha = 0.4;
%distance contour
surf(X_torus, Y_torus, Z_torus, 'FaceColor', 'r', 'FaceAlpha', falpha, 'edgecolor', 'none');
surf(X_sphere, Y_sphere, Z_sphere, 'FaceColor', 'r', 'FaceAlpha', falpha, 'edgecolor', 'none');
%

patch('XData', circ(1,:), 'YData', circ(2, :), 'ZData', z_circ, 'FaceColor', 'r', 'FaceAlpha', falpha, 'edgecolor', 'none');


%unsafe set
surf(rscale*X_sphere, rscale*Y_sphere, rscale*Z_sphere, 'FaceColor', 'r');
patch(circ(1,:), circ(2, :), zeros(size(phi)), 'FaceColor', 'r');
%unsafe set
% surf(X_sphere, Y_sphere, Z_sphere, 'FaceColor', 'r', 'FaceAlpha', falpha);
% patch(circ(1,:), circ(2, :), (R)*ones(size(phi)), 'FaceColor', 'r', 'FaceAlpha', falpha);

view(3);

daspect([1 1 1]) % preserves the shape of torus
% colormap('jet') % change color appearance 
title(['Distance Contour to half-sphere: ', num2str(r)], 'FontSize', 16)
xlabel('X');ylabel('Y');zlabel('Z');

