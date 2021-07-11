
%flow system with time-dependent uncertainty w (no structure)
SOLVE = 1;
PLOT = 0;

if SOLVE
    mset clear
%class-based flow distance implementation
mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('y', 2, 1)
mpol('w', 1, 1)

vars = struct;
vars.t = t;
vars.x = x;
vars.y = y;
vars.w = w;

lsupp = unsafe_support(vars);
BOX = 3;
lsupp = lsupp.set_box(BOX);
lsupp.Tmax = 5;

%initial set
C0 = [1.5; 0];
R0 = 0.4;
lsupp.X_init = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%unsafe set
% theta_c = pi/2;
theta_c =5*pi/4;
% Cu = [-0.5; -0.75];
Cu = [0; -0.7];
Ru = 0.5;
c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;

% theta_c = 3*pi/2;
% theta_c = 
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 
lsupp.X_unsafe = [c1f; c2f] >= 0;

lsupp.dist = (x-y)'*(x-y);

lsupp.disturb = (w^2 <= 1);
%% call distance manager

%flow system
wmax = 0.25;
draw = wmax*w;
% draw = 0;
f = [x(2); (-1 + draw)*x(1) + (1/3).* x(1).^3 - x(2) ];
PM = distance_manager(lsupp, f);

% loc = location_distance(lsupp, f);

order = 5;
d = 2*order;
% [objective, mom_con, supp_con] = PM.cons(d);
sol = PM.run(order);
dist_rec = sqrt(sol.obj_rec)
[optimal, mom_out, corner] = PM.loc.recover();
end

if PLOT
    load('flow_distance_noisy.mat')
    
    
    figure(53)
    clf
    hold on


    %distance contour
    x_dist_align = dist_contour(100, Ru, dist_rec);

    theta_cf = theta_c - 3*pi/2;
    Rot_mat = [cos(theta_cf) -sin(theta_cf); sin(theta_cf) cos(theta_cf)];
    x_dist = Rot_mat*x_dist_align + Cu;
    
    plot(x_dist(1, :), x_dist(2, :), 'r', 'DisplayName', 'Distance Contour', 'LineWidth', 2)
    
    
    theta = linspace(0, 2*pi, 200);
    circ = [cos(theta); sin(theta)];
    X0 = C0 + circ*R0;

    theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
    circ_half = [cos(theta_half_range); sin(theta_half_range)];
    Xu = Cu + circ_half* Ru;

%     out_sim = out_sim_multi;
    for i = 1:length(out_sim)
        if i == 1
            plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
        else
            plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'Off');
        end
    end

    plot(X0(1, :), X0(2, :), 'k', 'Linewidth', 3, 'DisplayName', 'Initial Set')
    patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')

    
    
       xlim([-1, 2.5])
    ylim([-1.5, 1.25])
%     axis square
    pbaspect([diff(xlim), diff(ylim), 1])
end

function x_dist = dist_contour(Ntheta, R, c)
    %compute a contour at distance c away from the half-circle with N_theta
    %sample points


    theta_q1 = linspace(0, pi/2, Ntheta);
    theta_q2 = linspace(pi/2, pi, Ntheta);
    theta_q34 = linspace(pi, 2*pi, 2*Ntheta);

    %contour level
    

    x_right = [c*cos(theta_q1)+R; c*sin(theta_q1)];
    x_left = [c*cos(theta_q2)-R; c*sin(theta_q2)];
    x_bottom = [(c+R)*cos(theta_q34); (c+R)*sin(theta_q34)];
    x_dist = [x_right, x_left, x_bottom];
end
% [objective, cons_eq, cons_ineq] = PM.loc.all_cons(d);