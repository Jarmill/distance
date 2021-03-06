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

%% call distance manager

%flow system
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
PM = distance_manager(lsupp, f);
% loc = location_distance(lsupp, f);

order = 4;
d = 2*order;
% [objective, mom_con, supp_con] = PM.cons(d);
[sol, PM] = PM.run(order);
sqrt(sol.obj_rec)
[optimal, mom_rec, corner_rec] = PM.loc.recover()
% [objective, cons_eq, cons_ineq] = PM.loc.all_cons(d);