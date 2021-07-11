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

%% unsafe set (moon)
%plot the moon
h_in = 0.4;
h_out = 1;
x_moon = moon_base(h_in, h_out);
 

%hugging the curve
moon_center = [0.4;-0.4];
moon_theta = -pi/10;
moon_scale = 0.8;

%same coordinates as half-circle example
% moon_center = [0;-0.7];
% moon_theta = -pi/4;
% moon_scale = 0.5

moon_rot = [cos(moon_theta), sin(-moon_theta); sin(moon_theta), cos(moon_theta)];
x_moon_move = moon_rot*x_moon*moon_scale + moon_center;


%statistics of the moon
c_in = [0;0.5*(1/h_in - h_in)];
r_in = 0.5*(1/h_in + h_in);

c_out = [0;0.5*(1/h_out - h_out)];
r_out = 0.5*(1/h_out + h_out);

c_in_scale = moon_rot*c_in*moon_scale + moon_center;
c_out_scale = moon_rot*c_out*moon_scale + moon_center;

r_in_scale = moon_scale*r_in;
r_out_scale = moon_scale*r_out;

%constraints of the moon
con_inner =  sum((y-c_in_scale).^2) - r_in_scale^2;
con_outer =  -sum((y-c_out_scale).^2) + r_out_scale^2;

lsupp.X_unsafe = [con_inner; con_outer] >= 0;

lsupp.dist = (x-y)'*(x-y);

%% call distance manager

%flow system
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
PM = distance_manager(lsupp, f);
% loc = location_distance(lsupp, f);

order = 5;
d = 2*order;
% [objective, mom_con, supp_con] = PM.cons(d);
sol = PM.run(order);
dist_rec = sqrt(sol.obj_rec)

[optimal, mom_rec, corner_rec] = PM.loc.recover();
% [objective, cons_eq, cons_ineq] = PM.loc.all_cons(d);

[recover, mom_rec, corner_rec] = PM.loc.recover();