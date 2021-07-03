%Perform distance estimation along flow system
%trajectories start from a circle X0, unsafe set is a half-circle Xu
%are trajectories safe with respect to the unsafe set?
%measure distance using the L1 norm
%Author: Jared Miller, 02/07/21
SAFE = 1;
order = 3; 
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


%angle of unsafe set 
if SAFE
%     theta_c = 5*pi/4;       
%   safe,   order 3: 0.3889 (suboptimal but still > 0)
%           objective: -3.8889488806487616e-01, mu =  1.8500490586470547e-24
    w_c = [-1; -1];
else
%     theta_c = 4*pi/4;       
%   unsafe, order 2: 1.2547e-08
%           objective: -1.5743918530608829e-16, mu = 6.6012067521220662e-39

    w_c = [-1; 0];
end
Cu = [0; -0.7];
%Cu = [2.5; 0];
Ru = 0.5; 

c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;
% w_c = [cos(theta_c); sin(theta_c)];
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

%L1 distance
%variables to implement distance
beta = sdpvar(2, 1);
consbeta = [beta >= -1; beta <= 1];

%distance penalty
c = sum(beta.*(x - y));


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
coeff = [cv; cw; gamma; coeff0; coeffaux; coeffdist; coeffL; beta];
cons = [cons0; consaux; consdist; consL; consbeta];

opts = sdpsettings('solver','mosek');
opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(cons, objective, opts, [coeff]);



dist_rec_l1 = value(gamma)