load('twist_4.mat')

dist = sqrt(sol.obj_rec);

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
figure(1)
clf

hold on


%unsafe set
surf(rscale*X_sphere + Cu(1), rscale*Y_sphere + Cu(2), rscale*Z_sphere+ Cu(3), 'FaceColor', 'r');
patch(circ(1,:)+ Cu(1), circ(2, :)+ Cu(2), zeros(size(phi)+ Cu(3)), 'FaceColor', 'r');

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

%% distance contour goes here


%Find the point-set L4 distance between a point x and a point y on a
%half-sphere in 3d. This problem is an SOCP
x = sdpvar(3, 1);

y = sdpvar(3, 1);

d = x-y; %distance to minimize

d2 = sdpvar(3, 1); %squared displacement
d4 = sdpvar(3, 1); %fourth power

Cu = [0.25; 0; 0];
Ru = 0.2;

c1f = (-(y(1)-Cu(1))^2 - (y(2)-Cu(2))^2 - (y(3)-Cu(3))^2 + Ru^2);
c2f = -y(3); 
cons = [c1f; c2f] >= 0;


M12 = cell(3, 1);
M24 = cell(3, 1);
for i = 1:3
    M12{i} = [1, d(i); d(i), d2(i)];
    M24{i} = [1, d2(i); d2(i), d4(i)];
    cons = [cons; M12{i} >= 0; M24{i} >= 0];
end

objective = sum(d4);

opts = sdpsettings('solver', 'mosek');
vars = [objective; y; d2; d4];
% sol = optimize(cons, objective, opts)
P = optimizer(cons, objective, opts, x, objective);
% sol = P(x0)
P_vec = @(x1, x2, x3) max(0, P([x1; x2; x3]));

% levellist = linspace(0, 1, 50);
% levellist =  0.041115989882192;

%L4^4 bound
levellist =  2.857873504325048e-06;

P_vec_bound = @(x1, x2, x3) P([x1; x2; x3]) - levellist;

% BOX = 1.5;
falpha = 0.4;
Fc = fimplicit3(@(x1,x2, x3) arrayfun(P_vec_bound, x1,x2, x3), [-0.25, 0.75, -0.5, 0.5, -0.5, 0.25], 'EdgeColor', 'None', 'FaceColor', 'r', 'FaceAlpha', falpha)
% Fc = fcontour3(@(x1,x2, x3) arrayfun(P_vec, x1,x2, x3), [-0.25, 0.75, -0.5, 0.5, -0.5, 0.25], 'LevelList', levellist)

view(-8.257407017987877, 17.199043276682510);
axis off

