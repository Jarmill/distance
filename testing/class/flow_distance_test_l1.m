clear
mset clear

%class-based flow distance implementation
mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('y', 2, 1)

vars = struct;
vars.t = t;
vars.x = x;
vars.y = y;

lsupp = unsafe_support(vars);
lsupp = lsupp.set_box(3);
lsupp.Tmax = 5;

%initial set
C0 = [1.5; 0];
R0 = 0.4;
lsupp.X_init = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%unsafe set
% theta_c = pi/2;
theta_c =5*pi/4; %ground truth

% theta_c = 5*pi/4 - 1/8;
 
% Cu = [-0.5; -0.75];
Cu = [0; -0.7];
Ru = 0.5;
c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;

w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 
lsupp.X_unsafe = [c1f; c2f] >= 0;

lsupp.dist = {[1;-1]*(x(1)-y(1)), [1;-1]*(x(2) - y(2))};

%% call distance manager

%flow system
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
PM = distance_manager(lsupp, f);
% loc = location_distance(lsupp, f);

order = 4;
d = 2*order;
% [objective, mom_con, supp_con] = PM.cons(d);
sol = PM.run(order);
sol.obj_rec
% [objective, cons_eq, cons_ineq] = PM.loc.all_cons(d);

%% Recovery
[opt_init, s_init, corner_init] = PM.loc.init.recover();
[opt_term, s_term, corner_term] = PM.loc.term.recover();
[opt_wass, s_wass, corner_wass] = PM.loc.wass{1}.recover();

opt_result = struct('xp', s_wass.x, 'y', s_wass.y, 'tp', s_term.t, 'x0', s_init.x, 'Cu', Cu,...
    'Ru', Ru, 'R0', R0, 'C0', C0);