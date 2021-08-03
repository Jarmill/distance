%Perform distance estimation along flow system
%trajectories start from a circle X0, unsafe set is a half-circle Xu
%are trajectories safe with respect to the unsafe set?
%Author: Jared Miller, 03/08/21
%
%Use correlative sparsity to decompose the constraint c(x,y)-w(x)>=0 over
%the set X times Y
SAFE = 1;
order = 4; 
d =2*order;

f_func = @(x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
Tmax = 5;
BOX = 3;


%% variables and support sets
t = sdpvar(1,1);
x = sdpvar(2,1);

gamma = sdpvar(1,1);
%initial set
C0 = [1.5; 0];
R0 = 0.4;
% 
% C0 = [1.5; 0];
% R0 = 0.4;


%angle of unsafe set 
if SAFE
    theta_c = 5*pi/4;       %safe, order 4: 0.2831, order 5: 0.2832
else
    theta_c = 4*pi/4;       %unsafe, order 4: 1.176e-3
end
Cu = [0; -0.7];
%Cu = [2.5; 0];
Ru = 0.5; 

c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 


%initial set
Xinit = struct('ineq', -(x(1)-C0(1))^2 - (x(2)-C0(2))^2 + R0^2, 'eq', []);

%state
X = struct('ineq', BOX^2 - x.^2, 'eq', []);

%unsafe
Y = struct('ineq', [c1f; c2f], 'eq', []);
    
%time
Tcons = struct('ineq', t*(1-t), 'eq', []);

All_state = struct('ineq', [Tcons.ineq; X.ineq], 'eq', []);
All_joint = struct('ineq', [X.ineq; Y.ineq], 'eq', []);


%polynomials
[v, cv] = polynomial([t; x], d);
[w, cw] = polynomial(x, d);

%L2 distance

c = sum((x-y).^2);

%test sos program 

%% formulate constraints
Lv = jacobian(v, t) + Tmax * jacobian(v,x)*f_func(x);

v0 = replace(v, t, 0);


[p0, cons0, coeff0] = constraint_psatz(v0 - gamma, Xinit, [x], d);

[paux, consaux, coeffaux] = constraint_psatz(w - v, All_state, [t;x], d);

[pL, consL, coeffL] = constraint_psatz(Lv, All_state, [t; x], d);

%objective

objective = -gamma;

%% the sparse distance constraint
% [pdist, consdist, coeffdist] = constraint_psatz(c - w, All_joint, [x; y], d);

%the cliques are I1 = (x1 x2 y1) and I2 = (x2 y1 y2)
%the constraint partition is J1 = (1 2) and J2 = (3 4)

[p0_1, cp0_1] = polynomial([x; y(1)], d);
[p0_2, cp0_2] = polynomial([x(2); y], d);

%deal with changing degrees in the constraint definitions later

%constraints on X
[px_1, cpx_1] = polynomial([x; y(1)], d-2);
[px_2, cpx_2] = polynomial([x; y(1)], d-2);

%constraints on Y
[py_1, cpy_1] = polynomial([x(2); y], d-2);
[py_2, cpy_2] = polynomial([x(2); y], d-2);


mult_dist = p0_1 + p0_2 + X.ineq'*[px_1; px_2] + Y.ineq'*[py_1; py_2];
p_dist = c - w - (mult_dist); % == 0 to enforce that c - w >= 0 on X times Y
[p_dist_coeff, p_dist_monom] = coefficients(p_dist, [x; y]);
coeffdist = [cp0_1; cp0_2; cpx_1; cpx_2; cpy_1; cpy_2];
consdist = [sos(p0_1); sos(p0_2); sos(px_1); sos(px_2); sos(py_1); sos(py_2); p_dist_coeff==0]; 


%% package up
coeff = [cv; cw; gamma; coeff0; coeffaux; coeffdist; coeffL];
cons = [cons0; consaux; consdist; consL];
opts = sdpsettings('solver', 'mosek', 'savesolverinput', 1);
opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);



dist_rec = sqrt(value(gamma))