%find the minimum distance between a shape carried along trajectories of the 'flow' system and a
%half-circle unsafe set

shape_color = [179, 0, 255]/255;

rng(343, 'twister');
status_feas = 1;

%options 
SOLVE_DIST = 1;
SOLVE_FEAS = 0;

SAMPLE = 1;

PLOT_FLOW = 1;
PLOT_DIST = 0;

n = 2;
order = 4;
d = 2*order;


%% problem parameters
f_func = @(x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
Tmax = 5;
BOX = 3;

%initial set
C0 = [1.5; 0];
R0 = 0.4;
% 
% C0 = [1.5; 0];
% R0 = 0.4;


%unsafe set
% theta_c = 3*pi/2;     %safe, bound = 0.18951
theta_c = 5*pi/4;      %safe, bound = 0.3184
% theta_c = 7*pi/4;     %unsafe, bound = 4.6647\times 10^{-4}

% Cu = [0; -0.5];
Cu = [0; -0.7];
%Cu = [2.5; 0];
Ru = 0.5;

%Shape Set
% shape_angle = 0;
% shape_angle = pi/6;
shape_angle = 5*pi/12;
shape_side = 0.1;


%plotting
FS_title = 16;

if SOLVE_DIST
mset clear
mset('yalmip',true);
% mset(sdpsettings('solver', 'mosek'));
% mset(sdpsettings('solver', 'mosek', 'mosek.MSK_DPAR_BASIS_TOL_S', 1e-9, ...
%                 'mosek.MSK_DPAR_BASIS_TOL_X', 1e-9, 'mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED', 1e-10, ...
%                 'mosek.MSK_DPAR_INTPNT_TOL_PATH', 1e-6));
%     
mset(sdpsettings('solver', 'mosek'))
mpol('x', 2, 1);
% f = Tmax * [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
f = Tmax*f_func(x);
%initial set
%C0 = [1.2; 0];

X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%half-circle set

% Cu = [0; -0.5];

c1f = Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;

% theta_c = 3*pi/2;
% theta_c = 
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 




%% set up measures
mpol('t0', 1, 1);
mpol('x0', 2, 1);
mu0 = meas([t0; x0]);

mpol('t_occ', 1, 1);
mpol('x_occ', 2, 1);
mu_occ = meas([t_occ; x_occ]);

mpol('tp', 1, 1);
mpol('xp', 2, 1);
mup = meas([tp; xp]);

%wasserstein
mpol('xw', 2, 1)
mpol('xu', 2, 1)
eta = meas([xw; xu]);

%shape measure
mpol('ss', 2, 1)
mpol('xs', 2, 1);
mus = meas([xs; ss]);

unsafe_cons = subs_vars([c1f; c2f], x, xu);

X0_con = subs_vars((x(1)-C0(1))^2 + (x(2)-C0(2))^2, x, x0) <= R0^2;
% X0_con = (x0 == C0);

%% support constraints
supp_con = [t0 == 0; t_occ*(1-t_occ)>=0; tp*(1-tp) >=0;
    x_occ.^2 <= BOX^2; xp.^2 <= BOX^2; xs.^2 <= BOX^2;
    ss.^2 <= shape_side^2;
    X0_con;
    xw.^2 <= BOX^2; 
    unsafe_cons >= 0;
    ];

%% moment constraints

%liouville constraint
y0 = mom(mmon([t0; x0], d));
yp = mom(mmon([tp; xp], d));



v  = mmon([t_occ; x_occ], d);
f_occ = subs_vars(f, x, x_occ);
Ay = mom(diff(v, t_occ) + diff(v, x_occ)*f_occ); 

    Liou = Ay + (y0 - yp);
    Liou_con = Liou == 0;

%marginal between peak and shape
ysx = mom(mmon(xs, d));

ypx = mom(mmon(xp, d));

ywx = mom(mmon(xw, d));

peak_shape = (ypx - ysx == 0);

%marginal between shape and wass
R_shape = [cos(shape_angle) -sin(shape_angle); sin(shape_angle) cos(shape_angle)];
shape_transform = xs + R_shape*ss;
push_s = subs(mmon(xs, d), xs, shape_transform);

wass_shape = (ywx - mom(push_s) == 0);
    
%constraints together
mom_con = [Liou_con; wass_shape; peak_shape; mass(mu0)==1];


%objective
mom_dist = mom(sum((xw-xu).^2));
objective = min(mom_dist);

%% solve problem
%Input LMI moment problem
P = msdp(objective, ...
    mom_con, supp_con);

%solve LMIP moment problem
[status, obj, m, dual_rec] = msol(P);    

%% analyze solutions

dist_rec = sqrt(double(mom_dist));
disp(['distance bound: ', num2str(dist_rec)])
M0 = double(mmat(mu0));
% Mocc = double(mmat(mu_occ));
Mp = double(mmat(mup));
Mw = double(mmat(eta));
Ms = double(mmat(mus));

M0_1 = M0(1:(n+2), 1:(n+2));
Mp_1 = Mp(1:(n+2), 1:(n+2));
Mw_1 = Mw(1:(2*n+1), 1:(2*n+1));
Ms_1 = Ms(1:(2*n+1), 1:(2*n+1));

rankp = rank(Mp_1, 1e-3);
rank0 = rank(M0_1, 1e-3);
rankw = rank(Mw_1, 1e-3);
ranks = rank(Ms_1, 1e-3);

xu_rec = double(mom(xu));
xp_rec = double(mom(xw));
x0_rec = double(mom(x0));
xs_rec = double(mom(xs));
ss_rec = double(mom(ss));
tp_rec = Tmax*double(mom(tp));


optimal_pt = all([rankp; rank0; rankw; ranks]==1);

end 

%% Feasibility program (should probably be solved first)
if SOLVE_FEAS
    
    %initial in X0
    %terminal (peak) in Xu
    u_cons_feas = subs_vars([c1f; c2f], x, xp);
    
    supp_con_feas = [t0 == 0; t_occ*(1-t_occ)>=0; tp*(1-tp) >=0;
    x_occ.^2 <= BOX^2;
    X0_con;
    u_cons_feas >= 0];
%     xw.^2 <= BOX^2; 
%     u_cons >= 0;
    
    
    mom_con_feas = [Liou_con; mass(mu0)==1];
    objective_feas = min(mom(xp(2)));
    
    P_feas = msdp(objective_feas, ...
    mom_con_feas, supp_con_feas);

    %solve LMIP moment problem
    [status_feas, obj_feas, m_feas, dual_rec_feas] = msol(P_feas);    

%     M0_1_feas = M0(1:(n+2), 1:(n+2));
%     Mp_1_feas = Mp(1:(n+2), 1:(n+2));
% 
%     rankp_feas = rank(Mp_1, 1e-3);
%     rank0_feas = rank(M0_1, 1e-3);
% 
%     xp_feas_rec = double(mom(xp));
%     x0_feas_rec = double(mom(x0));
%     tp_feas_rec = Tmax*double(mom(tp));

end

%% Sample trajectories
if SAMPLE
    Nsample = 150;
    Tmax_sim = 5;
%     sampler = @() circle_sample(1)'*R0 + C0;

    
    flow_event = @(t, x) box_event(t, x, BOX);
    sample_x = @() circle_sample(1)'*R0 + C0;    

    
    ode_options = odeset('Events',flow_event, 'MaxStep', 0.05, 'AbsTol', 1e-7, 'RelTol', 1e-6);;

    out_sim = cell(Nsample, 1);
    
    %distance function
    dist_func = @(x_in) aff_half_circ_dist(x_in, Ru, theta_c, Cu);
%     peak_traj_dist = arrayfun(@(i) dist_func(out_sim_peak.x(i, :)'),...
%         1:size(out_sim_peak.x, 1));
    
    for i = 1:Nsample
        x0_curr = sample_x();        
        [time_curr, x_curr] = ode15s(@(t, x) f_func(x), [0, Tmax], x0_curr, ode_options);
        dist_curr = arrayfun(@(i) dist_func(x_curr(i, :)'),...
            1:size(x_curr, 1))';
        
        out_sim{i} = struct('t', time_curr, 'x', x_curr, 'dist', dist_curr);
        
        %figure out coordinate transformation to evaluate true distances
    end
    
    if optimal_pt
        [time_opt_traj, x_opt_traj] = ode15s(@(t, x) f_func(x), [0, Tmax], x0_rec, ode_options);
        dist_opt_traj = arrayfun(@(i) dist_func(x_opt_traj(i, :)'),...
            1:size(x_opt_traj, 1));
        out_sim_peak = struct('t',   time_opt_traj, 'x', x_opt_traj, 'dist', dist_opt_traj);
    end
    
%     s_opt = sampler_options;
%     s_opt.sample.x = @() circle_sample(1)'*R0 + C0;    
%     s_opt.Tmax = Tmax;

%     s_opt.parallel = 0;
%     s_opt.mu = 0.4;
%     out_sim = sampler(dynamics, Nsample, s_opt);
%     out.optimal = 0;

end

%% Plot flow output
if PLOT_FLOW
    
    %initial and unsafe sets
    theta = linspace(0, 2*pi, 100);
    circ = [cos(theta); sin(theta)];
    %half_theta = linspace(pi/2, 3*pi/2, 100);
    %half_circ = [cos(half_theta); sin(half_theta)];
    
    %initial set    
    X0 = C0 + circ*R0;
    
    %sample from X0 
    rng(33, 'twister')

            
    if theta_c > pi
        theta_c = theta_c - 2 * pi;
    end
    
    theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
    circ_half = [cos(theta_half_range); sin(theta_half_range)];
    Xu = Cu + circ_half* Ru;
    
    figure(1)
    clf
    hold on
    
    
%     for i = 1:Nsample
%         if i == 1
%             plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
%         else
%             plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'Off');
%         end
%     end
%     
    rect_shape_0 = R_shape*shape_side*[1,1,-1,-1,1;-1,1,1,-1,-1];
    rect_shape_init = rect_shape_0 + x0_rec;
    rect_shape_peak = rect_shape_0+ xs_rec;

    plot(X0(1, :), X0(2, :), 'k', 'Linewidth', 3, 'DisplayName', 'Initial Set')
    patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')
%distance contour
    x_dist_align = dist_contour(100, Ru, dist_rec);

    theta_cf = theta_c - 3*pi/2;
    Rot_mat = [cos(theta_cf) -sin(theta_cf); sin(theta_cf) cos(theta_cf)];
    x_dist = Rot_mat*x_dist_align + Cu;
    
%     plot(x_dist(1, :), x_dist(2, :), 'r', 'DisplayName', 'Distance Contour', 'LineWidth', 2)
    
    
    
    if optimal_pt

        plot(out_sim_peak.x(:, 1), out_sim_peak.x(:, 2), 'b', 'DisplayName', 'Closest Traj.', 'LineWidth', 2);       
        
        scatter(x0_rec(1), x0_rec(2), 200, 'ob', 'DisplayName', 'Closest Initial', 'LineWidth', 2);        
        scatter(xs_rec(1), xs_rec(2), 600, '.b', 'DisplayName', 'Closest Shape Center', 'LineWidth', 2);        
        
        scatter(xp_rec(1), xp_rec(2), 200, '*b', 'DisplayName', 'Closest Point', 'LineWidth', 2);        
        scatter(xu_rec(1), xu_rec(2), 200, 'sb', 'DisplayName', 'Closest Unsafe', 'LineWidth', 2);        
        
         patch(rect_shape_init(1, :), rect_shape_init(2, :),'k', 'Linewidth', 3,'DisplayName', 'Initial Shape', 'EdgeColor', shape_color, 'FaceColor', 'None')
    
        patch(rect_shape_peak(1, :), rect_shape_peak(2, :),'k', 'Linewidth', 3,'DisplayName', 'Peak Shape', 'EdgeColor', shape_color, 'FaceColor', 'None')
        plot([xp_rec(1); xu_rec(1)], [xp_rec(2); xu_rec(2)], ':k', 'DisplayName', 'Closest Distance', 'Linewidth', 1.5)
       
    end

    
    legend('location', 'northwest')
    
%     xlim([-0.6, 1.7])
%     ylim([-1.3, 0.3])
    xlim([-0.6, 1.9])
    ylim([-1.3, 0.5])
    pbaspect([diff(xlim), diff(ylim), 1])
    xlabel('x_1')
    ylabel('x_2')
%     axis square
    
    if status_feas == 0
        %feasible program 
        title_str = (['Unsafe, false L_2 distance bound is ', num2str(dist_rec, 3)]);
    else    
        title_str = (['L_2 distance bound is ', num2str(dist_rec, 3)]);        
    end
    title(title_str, 'FontSize' , FS_title)
end

%% Distance plot
if PLOT_DIST
    
    
    
%     peak_traj_dist = arrayfun(@(i) half_circ_dist(out_sim_peak.x(i, :)'-Cu, Ru),...
%         1:size(out_sim_peak.x, 1));
%     
    
    figure(2)
    clf
    hold on
    for i = 1:Nsample
        if i == 1
            plot(out_sim{i}.t, out_sim{i}.dist, 'c', 'DisplayName', 'Trajectories');
        else
            plot(out_sim{i}.t, out_sim{i}.dist, 'c', 'HandleVisibility', 'Off');
        end
    end
    
    plot(xlim, dist_rec*[1,1], '--r', 'LineWidth', 2 ,'DisplayName', 'Distance Bound');
    
    if optimal_pt
        plot(out_sim_peak.t, out_sim_peak.dist, 'b', 'LineWidth', 2, 'DisplayName','Closest Traj.');
        scatter(tp_rec, dist_rec, 300, '*b', 'LineWidth', 2,'DisplayName','Closest Point');
    end
    
    xlabel('time')
    ylabel('distance to unsafe set')
    title('Distance to unsafe set along trajectories', 'FontSize' , FS_title)
    legend('location', 'northwest')
end


function dist_out = half_circ_dist(x_in, R)
    %return the L2 distance  between the point x_in and the half circle
    %||x_in||^2 <= R^2 intersect x_in(2) <= 0.
%     reshape(x_in, [], 1);
    if x_in(2) >= 0
        %flat region
        if x_in(1) < -R
            dist_out = hypot(x_in(1)+R, x_in(2));
        elseif x_in(1) > R
            dist_out = hypot(x_in(1)-R, x_in(2));
        else
            dist_out = x_in(2);
        end
    else
        %circle region
        dist_out = max(norm(x_in, 2)-R, 0);
    end

end

function dist_out = aff_half_circ_dist(x_in, R, theta_c, Cu)

    theta_cf = theta_c - 3*pi/2;
    Rot_mat = [cos(theta_cf) -sin(theta_cf); sin(theta_cf) cos(theta_cf)];
    x_aff = Rot_mat'*(x_in - Cu);
    
    dist_out = half_circ_dist(x_aff, R);

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
