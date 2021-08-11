%uses symmetric routines from hybrid_peak
%% properties of system
BOX = 3;
shape_size = 0.1;
%Tmax = 5;
Tmax = 3;

order = 1;
d = 2*order;
% ang_vel = 1;
ang_vel = 0;
%dynamics
% f_func = @(x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
f_func = @(om) [-ang_vel*om(2); ang_vel*om(1); om(3); -om(3) + (1/3).* om(3).^3 - om(4) ];

%perform distance estimation 

%% variable definiftion
%sym order 2: 
%std order 2:   c
%sym:
%std:           s, x, y, sx, sy
t = sdpvar(1, 1);  %time 
om = sdpvar(4, 1); %orientation (cos(theta), sin(theta), x, y)
x = sdpvar(2, 1);  %point in state
y = sdpvar(2, 1);  %point on unsafe set
s = sdpvar(2, 1);  %point on shape

gamma = sdpvar(1, 1);

% var_count = [0; 1; 0; 5];
% mom_out = sym_moments(var_count, 2);

%% support sets

%initial set
C0 = [1.5; 0];
R0 = 0.4;
Xinit = struct('ineq', -(om(3)-C0(1))^2 - (om(4)-C0(2))^2 + R0^2, 'eq', []);
theta0 = [cos(5*pi/12); sin(5*pi/12)];
om0 = [theta0; [1.48887368624732;-0.399810603537798]];

%state 


%orientation set
Om = struct('ineq', BOX^2 - om(3:4).^2, 'eq', [om(1)^2 + om(2)^2 - 1]);

%total set
SOm = struct('ineq', [Om.ineq; 1 - s.^2], 'eq', Om.eq);

%unsafe set
theta_c = 5*pi/4;       %safe, order 4: 0.2831, order 5: 0.2832
Cu = [0; -0.7];
Ru = 0.5; 
c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 

Y = struct('ineq', [c1f; c2f], 'eq', []);
X = struct('ineq', BOX^2 - x.^2, 'eq', []);

XY = struct('ineq', [X.ineq; Y.ineq], 'eq', [X.eq; Y.eq]);

%time set
Tcons = struct('ineq', t*(1-t), 'eq', []);

TOm = struct('ineq', [Om.ineq; Tcons.ineq], 'eq', Om.eq);

%% polynomials
%v(t, omega)
var_count_v = [0,1,0,4];
mom_out_v = sym_moments(var_count_v, order);
v_monom = recovermonoms(mom_out_v.monom_int,[om; t]);
cv = sdpvar(length(v_monom), 1);
v = cv'*v_monom;

%z(omega)
var_count_z = [0,1,0,3];
mom_out_z = sym_moments(var_count_z, order);
z_monom = recovermonoms(mom_out_z.monom_int,om);
cz = sdpvar(length(z_monom), 1);
z = cz'*z_monom;

%w(x)
[w, cw] = polynomial(x, d);


%% Initial constraint
% v0 = replace(v, [om(1:2);t], [theta0;0]);
% [p_init, cons_init, Gram_init] = sym_yalmip_psatz(v0 - gamma, Xinit, order, om(3:4), [0,0,0,2]);

v0 = replace(v, [om; t], [om0; 0]);
cons_init = (v0 - gamma >= 0);

%% Lie Derivative
Lv = Tmax*jacobian(v, om)*f_func(om) + jacobian(v, t);
[p_lie, cons_lie, Gram_lie] = sym_yalmip_psatz(Lv, TOm, order, [om; t], [0,1,0,4]);


%% orientation and proxy (z and v)
[p_zv, cons_zv, Gram_zv] = sym_yalmip_psatz(z-v, TOm, order, [om; t], [0,1,0,4]);

%% shape relation (w and z)
%this is likely the bottleneck
shape_transform = om(3:4) + [om(1), -om(2); om(2), om(1)]*(shape_size*s);
shape_transform = om(3:4);
w_push = replace(w, x, shape_transform);

% [p_wz, cons_wz, Gram_wz] = sym_yalmip_psatz(w_push - z, SOm, 2*order, [om; s], [0,1,0,5]);
[p_wz, cons_wz, Gram_wz] = sym_yalmip_psatz(w_push - z, Om, order, [om], [0,1,0,3]);

%% cost proxy
%this can be decomposed through csp if desired
c = sum((x-y).^2);
[p_dist, cons_dist, Gram_dist] = sym_yalmip_psatz(c-w, XY, order, [x; y], [0,0,0,4]);


cons = [cons_init:'init'; cons_lie:'lie'; cons_zv:'zv'; cons_wz:'wz'; cons_dist:'cw'];
objective = -gamma;



opts = sdpsettings('solver', 'mosek', 'savesolverinput', true);
sol = optimize(cons, objective, opts);

% value(gamma)
dist_rec = sqrt(value(gamma));
fprintf('rotating shape L2 bound = %0.4f\n', dist_rec)
