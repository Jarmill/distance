%Certify set disconnectedness
%union of [0, low] and [high, 1]

order = 5; 
d =2*order;

f_func = @(x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
Tmax = 5;
BOX = 3;


%% variables and support sets
t = sdpvar(1,1);
x = sdpvar(2,1);
y = sdpvar(2,1);

gamma = sdpvar(1,1);
%initial set
C0 = [1.5; 0];
R0 = 0.4;
% 
% C0 = [1.5; 0];
% R0 = 0.4;


%unsafe set
% theta_c = 5*pi/4;       %safe, order 4: 0.2831, order 5: 0.2832
theta_c = 4*pi/4;     %unsafe, bound = 4.6647\times 10^{-4}

Cu = [0; -0.7];
%Cu = [2.5; 0];
Ru = 0.5; 

c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;

% theta_c = 3*pi/2;
% theta_c = 
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

[pdist, consdist, coeffdist] = constraint_psatz(c - w, All_joint, [x; y], d);

[pL, consL, coeffL] = constraint_psatz(Lv, All_state, [t; x], d);

%objective

objective = -gamma;

%% package up
coeff = [cv; cw; gamma; coeff0; coeffaux; coeffdist; coeffL];
cons = [cons0; consaux; consdist; consL];
opts = sdpsettings('solver', 'mosek');
opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);

sqrt(value(gamma))