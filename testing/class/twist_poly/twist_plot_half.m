load('twist_5.mat')

% dist = sqrt(sol.obj_rec);

%% compute the initial set
[Xs, Ys, Zs] = sphere(40);
X0 = R0*Xs + C0(1);
Y0 = R0*Ys + C0(2);
Z0 = R0*Zs + C0(3);

%% compute the unsafe set

%quarter-torus
R = Ru;  %radius of half-sphere
r = dist; %inner radius of torus

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


%% plot the sets
figure(2)
clf

hold on
falpha = 0.4;
%distance contour
surf(X_torus+ Cu(1), Y_torus+ Cu(2), Z_torus+ Cu(3), 'FaceColor', 'r', 'FaceAlpha', falpha, 'edgecolor', 'none');
surf(X_sphere+ Cu(1), Y_sphere+ Cu(2), Z_sphere+ Cu(3), 'FaceColor', 'r', 'FaceAlpha', falpha, 'edgecolor', 'none');
%

patch('XData', circ(1,:)+ Cu(1), 'YData', circ(2, :)+ Cu(2), 'ZData', z_circ+ Cu(3), 'FaceColor', 'r', 'FaceAlpha', falpha, 'edgecolor', 'none');


%unsafe set
surf(rscale*X_sphere + Cu(1), rscale*Y_sphere + Cu(2), rscale*Z_sphere+ Cu(3), 'FaceColor', 'r');
patch(circ(1,:)+ Cu(1), circ(2, :)+ Cu(2), zeros(size(phi)+ Cu(3)), 'FaceColor', 'r');
%unsafe set
% surf(X_sphere, Y_sphere, Z_sphere, 'FaceColor', 'r', 'FaceAlpha', falpha);
% patch(circ(1,:), circ(2, :), (R)*ones(size(phi)), 'FaceColor', 'r', 'FaceAlpha', falpha);

view(3);


% colormap('jet') % change color appearance 
% title(['Distance Contour to half-sphere: ', num2str(r)], 'FontSize', 16)
% xlabel('X');ylabel('Y');zlabel('Z');

%% plot trajectories
surf(X0, Y0, Z0, 'FaceColor', 0.5*[1,1,1], 'FaceAlpha', falpha, 'edgecolor', 'none');

for i = 1:length(out_sim)
    xt = out_sim{i}.x;
plot3(xt(:, 1), xt(:, 2), xt(:, 3), 'c')
% scatter3(xt(1, 1), xt(1, 2), xt(1, 3), 100, 'k')
end

daspect([1 1 1]) % preserves the shape of torus
view(-8.257407017987877, 17.199043276682510);