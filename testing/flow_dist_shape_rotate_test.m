%find the minimum distance between a shape carried along trajectories of the 'flow' system and a
%half-circle unsafe set

%the shape now rotates along trajectories 
%at order 2, the shape moment matrix has size 210
%at order 3, the shape moment matrix has size 924

% quotient structure in w(4)^2 = 1- w(3)^2 and for all measures including
% w. Also term sparsity everywhere, but most particularly in the shape
% measure. The shape measure for this example has moments up to degree 4d.

shape_color = [179, 0, 255]/255;

rng(343, 'twister');
status_feas = 1;

%options 
SOLVE_DIST = 0;
SOLVE_FEAS = 0;

SAMPLE = 0;

PLOT_FLOW = 1;
PLOT_DIST = 0;

MOM_SUB = 1;

n = 2;
order =3;
d = 2*order;


%% problem parameters
ang_vel = 1;
%w [horizontal coord; vertical coord; cosine angle; sine angle]
f_func = @(w) [w(2); -w(1) + (1/3).* w(1).^3 - w(2); -w(4)*ang_vel; w(3)*ang_vel];
Tmax = 5;
BOX = 3;

%initial set
C0 = [1.5; 0];
R0 = 0.4;


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
mpol('w', 4, 1);
% f = Tmax * [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
f = Tmax*f_func(w);
%initial set
%C0 = [1.2; 0];

% X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);
% W0 = [(w(1)-C0(1))^2 + (w(2)-C0(2))^2 <= R0^2); w(3) ==cos(shape_angle); w(4) == sin(shape_angle)];
%half-circle set

% Cu = [0; -0.5];

c1f = Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;

% theta_c = 3*pi/2;
% theta_c = 
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 




%% set up measures
mpol('t0', 1, 1);
mpol('w0', 4, 1);
mu0 = meas([t0; w0]);

mpol('t_occ', 1, 1);
mpol('w_occ', 4, 1);
mu_occ = meas([t_occ; w_occ]);

mpol('tp', 1, 1);
mpol('wp', 4, 1);
mup = meas([tp; wp]);

%wasserstein
mpol('xx', 2, 1)
mpol('xu', 2, 1)
eta = meas([xx; xu]);

%shape measure
mpol('ss', 2, 1)
mpol('ws', 4, 1);
mus = meas([ws; ss]);

unsafe_cons = subs_vars([c1f; c2f], x, xu);

% X0_con = subs_vars((x(1)-C0(1))^2 + (x(2)-C0(2))^2, x, x0) <= R0^2;
% X0_con = subs_vars((x(1)-C0(1))^2 + (x(2)-C0(2))^2, x, x0) <= R0^2;
% X0_con = (x0 == C0);
W0_con = [((w0(1)-C0(1))^2 + (w0(2)-C0(2))^2 <= R0^2); (w0(3)==cos(shape_angle)); (w0(4) == sin(shape_angle))];
% circ_con = [w(3)^2 + w(4)^2 == 1];

Wall_con = [w(1:2).^2 <= BOX^2; w(3)^2 + w(4)^2 == 1];


%% support constraints
supp_con = [t0 == 0; t_occ*(1-t_occ)>=0; tp*(1-tp) >=0;    
    subs(Wall_con, w, w_occ);
    subs(Wall_con, w, wp);
    subs(Wall_con, w, ws);
    ss.^2 <= 1;
    W0_con;
    xx.^2 <= BOX^2; 
    unsafe_cons >= 0; %xu
    ];

%% moment constraints

%liouville constraint
y0 = mom(mmon([t0; w0], d));
yp = mom(mmon([tp; wp], d));



v  = mmon([t_occ; w_occ], d);
f_occ = subs_vars(f, w, w_occ);
Ay = mom(diff(v, t_occ) + diff(v, w_occ)*f_occ); 

if MOM_SUB
    Liou_con = (yp == Ay + y0);
else
    Liou = Ay + (y0 - yp);
    Liou_con = (Liou == 0);
end
%marginal between peak and shape
ysw = mom(mmon(ws, d));

ypw = mom(mmon(wp, d));

ywx = mom(mmon(xx, d));

%

if MOM_SUB
     peak_shape = (ysw == ypw);
else
     peak_shape = (ypw - ysw == 0);
end

%marginal between shape and wass
% R_shape = [cos(shape_angle) -sin(shape_angle); sin(shape_angle) cos(shape_angle)];
% shape_transform = xs + R_shape*ss;
shape_transform = ws(1:2) + [ws(3), -ws(4); ws(4), ws(3)]*ss*shape_side;
push_s = subs(mmon(xx, d), xx, shape_transform);

% if MOM_SUB
%     wass_shape = (ywx == mom(push_s));
% else
    wass_shape = (ywx - mom(push_s) == 0);
% end
%constraints together
mom_con = [Liou_con; wass_shape; peak_shape; mass(mu0)==1];


%objective
mom_dist = mom(sum((xx-xu).^2));
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

M0_1 = double(mom([1;w0]*[1;w0]'));
Mp_1 = double(mom([1;tp;wp]*[1;tp;wp]'));
Mw_1 = double(mom([1;xx;xu]*[1;xx;xu]'));
Ms_1 = double(mom([1;ss;ws]*[1;ss;ws]'));



rankp = rank(Mp_1, 1e-3);
rank0 = rank(M0_1, 1e-3);
rankw = rank(Mw_1, 1e-3);
ranks = rank(Ms_1, 1e-3);

xu_rec = double(mom(xu));
xp_rec = double(mom(xx));
w0_rec = double(mom(w0));
ws_rec = double(mom(ws));
ss_rec = double(mom(ss));
tp_rec = Tmax*double(mom(tp));


optimal_pt = all([rankp; rank0; rankw; ranks]==1);


rot_shape = struct('xu', xu_rec, 'xp', xp_rec, 'w0', w0_rec, 'ws', ws_rec,...
    'ss', ss_rec, 'tp', tp_rec, 'M0', M0_1, 'Mp', Mp_1, 'Ms', Ms_1, 'Mw', Mw_1,...
    'order', order, 'dist', dist_rec);


end 


%% Sample trajectories
if SAMPLE
    Nsample = 150;
    Tmax_sim = 5;
%     sampler = @() circle_sample(1)'*R0 + C0;

    
    flow_event = @(t, x) box_event(t, x, BOX);
    sample_x = @() [circle_sample(1)'*R0 + C0;  [cos(shape_angle); sin(shape_angle)];];  

    
    ode_options = odeset('Events',flow_event, 'MaxStep', 0.05, 'AbsTol', 1e-7, 'RelTol', 1e-6);;

    out_sim = cell(Nsample, 1);
    
    %distance function
%     dist_func = @(x_in) aff_half_circ_dist(x_in, Ru, theta_c, Cu);
%     peak_traj_dist = arrayfun(@(i) dist_func(out_sim_peak.x(i, :)'),...
%         1:size(out_sim_peak.x, 1));
    
    for i = 1:Nsample
        x0_curr = sample_x();        
        [time_curr, x_curr] = ode15s(@(t, x) f_func(x), [0, Tmax], x0_curr, ode_options);
%         dist_curr = arrayfun(@(i) dist_func(x_curr(i, :)'),...
%             1:size(x_curr, 1))';
        dist_curr = [];
        
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

    R_shape_0 = [w0_rec(3), -w0_rec(4); w0_rec(4), w0_rec(3)];
    rect_shape_0 = R_shape_0*shape_side*[1,1,-1,-1,1;-1,1,1,-1,-1];
    
    R_shape_p = [ws_rec(3), -ws_rec(4); ws_rec(4), ws_rec(3)];
    rect_shape_p = R_shape_p*shape_side*[1,1,-1,-1,1;-1,1,1,-1,-1];
    rect_shape_init = rect_shape_0 + w0_rec(1:2);
    rect_shape_peak = rect_shape_p+ ws_rec(1:2);

    plot(X0(1, :), X0(2, :), 'k', 'Linewidth', 3, 'DisplayName', 'Initial Set')
    patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')
%distance contour
    x_dist_align = dist_contour(100, Ru, dist_rec);

    theta_cf = theta_c - 3*pi/2;
    Rot_mat = [cos(theta_cf) -sin(theta_cf); sin(theta_cf) cos(theta_cf)];
    x_dist = Rot_mat*x_dist_align + Cu;
    
%     plot(x_dist(1, :), x_dist(2, :), 'r', 'DisplayName', 'Distance Contour', 'LineWidth', 2)
    
    
    
    if optimal_pt
        flow_event = @(t, x) box_event(t, x, BOX);
            ode_options = odeset('Events',flow_event, 'MaxStep', 0.05, 'AbsTol', 1e-7, 'RelTol', 1e-6);;

        [time_opt_traj, w_opt_traj] = ode15s(@(t, w) f_func(w), ...
            [0, Tmax], w0_rec, ode_options);

        plot(w_opt_traj(:, 1), w_opt_traj(:, 2), 'b', 'DisplayName', 'Closest Traj.', 'LineWidth', 2);       
        
        scatter(w0_rec(1), w0_rec(2), 200, 'ob', 'DisplayName', 'Closest Initial', 'LineWidth', 2);        
        scatter(ws_rec(1), ws_rec(2), 600, '.b', 'DisplayName', 'Closest Shape Center', 'LineWidth', 2);        
        

        patch(rect_shape_init(1, :), rect_shape_init(2, :),'k', 'Linewidth', 3,'DisplayName', 'Initial Shape', 'EdgeColor', shape_color, 'FaceColor', 'None')
    
        patch(rect_shape_peak(1, :), rect_shape_peak(2, :),'k', 'Linewidth', 3,'DisplayName', 'Peak Shape', 'EdgeColor', shape_color, 'FaceColor', 'None')
        plot([xp_rec(1); xu_rec(1)], [xp_rec(2); xu_rec(2)], ':k', 'DisplayName', 'Closest Distance', 'Linewidth', 1.5)
       
        scatter(xp_rec(1), xp_rec(2), 200, '*b', 'DisplayName', 'Closest Point', 'LineWidth', 2);        
        scatter(xu_rec(1), xu_rec(2), 200, 'sb', 'DisplayName', 'Closest Unsafe', 'LineWidth', 2);        
        
        
    end

    
%     legend('location', 'northwest')
    
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
    
    
    title('')
    axis off
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
