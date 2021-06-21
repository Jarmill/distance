%% variable definition
SOLVE = 0;
PLOT = 1;

if SOLVE
%state (time-independent)
x = sdpvar(2,1);

%dynamics
f = [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];

%% support sets
%total set
BOX = 3;
X = struct('ineq', BOX - x.^2, 'eq', []);

%initial set

%initial set
C0 = [1.5; 0];
R0 = 0.4;
X_init = struct('ineq', R0^2 - ((x(1)-C0(1))^2 + (x(2)-C0(2))^2), 'eq', []);

%unsafe set
theta_c =5*pi/4; %ground truth
Cu = [0; -0.7];
Ru = 0.5;
c1f = Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;

w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 
X_unsafe = struct('ineq', [c1f; c2f], 'eq', 0);

%% Polynomials

order = 3;
d = 2*order;

[B, cB] = polynomial(x, d);

LieB = jacobian(B, x)*f;

%B >= epsilon on X_init
%B < 0 on X_unsafe
epsilon = 0.01;

%% Constraints

% [p, cons, coeff] = constraint_psatz(f, X, vars, d);        [p, cons, coeff] = constraint_psatz(f, X, vars, d);        
[p0, con0, coeff0] = constraint_psatz(B-epsilon, X_init, x, d);        
[pu, conu, coeffu] = constraint_psatz(-B, X_unsafe, x, d);        

[pf, conf, coefff] = constraint_psatz(LieB, X, x, d);        


con = [con0; conu; conf];
coeff = [coeff0; coeffu; coefff; cB];

%% solve barrier system
opts = sdpsettings('solver', 'mosek', 'verbose', 1);
opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(con, 0, opts, coeff);


%% recover solution
[cBr, mBr] = coefficients(B, x);
B_rec = value(cBr)'*mBr;

B_func = polyval_func(B_rec, x);
end

%% plot solution
% out_sim = load('flow_distance.mat', 'out_sim');
if PLOT
load('flow_distance.mat', 'out_sim');
figure(30)
clf
hold on 
for i = 1:length(out_sim)
    if i == 1
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
    else
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'Off');
    end
end
plot(X0(1, :), X0(2, :), 'k', 'Linewidth', 3, 'DisplayName', 'Initial Set')
patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')

fimplicit(@(x, y) B_func([x;y]), '--', 'Color', [ 0.4416    0.7490    0.4322], 'LineWidth', 3)

       xlim([-1, 2.5])
    ylim([-1.5, 1.25])
%     axis square
    pbaspect([diff(xlim), diff(ylim), 1])
    axis off
end