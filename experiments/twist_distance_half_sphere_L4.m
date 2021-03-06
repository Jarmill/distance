clear
mset clear

%class-based flow distance implementation
mpol('t', 1, 1)
mpol('x', 3, 1)
mpol('y', 3, 1)

vars = struct;
vars.t = t;
vars.x = x;
vars.y = y;

lsupp = unsafe_support(vars);
lsupp = lsupp.set_box(1);
lsupp.Tmax = 5;

%initial set
C0 = [-0.5; 0; 0];
R0 = 0.2;
lsupp.X_init = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2  + (x(3) - C0(3))^2 <= R0^2);

%unsafe set
Cu = [0.25; 0; 0];
Ru = 0.2;

c1f = (-(y(1)-Cu(1))^2 - (y(2)-Cu(2))^2  - (y(3) - Cu(3))^2 + Ru^2);

c2f = -y(3); 
lsupp.X_unsafe = [c1f; c2f] >= 0;

%L2 bound
% lsupp.dist = (x-y)'*(x-y);
% lsupp.dist = sum((x-y).^2);

lsupp.dist = sum((x-y).^4);

lsupp.CSP = 1;

%% call distance manager

%flow system
% f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
%twist dynamics
A_true = [-1 1 1; -1 0 -1; 0 1 -2];
B_true = [-1 0 -1;
          0 1 1;
          1 1 0]/2;

f_true = @(t,x) A_true*x - B_true*(4*x.^3 - 3*x);
f = f_true(0, x);
PM = distance_manager(lsupp, f);
% loc = location_distance(lsupp, f);

order = 1;
d = 2*order;
% [objective, mom_con, supp_con] = PM.cons(d);
sol = PM.run(order);
fprintf('L4 bound: %0.5f\n', (sol.obj_rec)^(0.25))
% fprintf('L2 bound: %0.5f\n', sqrt(sol.obj_rec))
% [objective, cons_eq, cons_ineq] = PM.loc.all_cons(d);

[recover, mom_rec, corner_rec] = PM.loc.recover();