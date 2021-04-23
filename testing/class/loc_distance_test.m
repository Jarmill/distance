mset clear
%% variables
BOX = 1;

mpol('t', 1, 1)
mpol('x', 2, 1)
mpol('th', 1, 1)
mpol('w', 2, 1)
mpol('y', 2, 1)


vars = struct;
vars.t = t;
vars.x = x;
vars.th = th;
vars.w = w;
vars.y = y;


if BOX
    mpol('b', 2, 1)
    vars.b = b;
end

%% location support 
X = x.^2 <= 1;
TH = th.^2 <= 1;
W = sum(w.^2) <= 1;

X_init = (x-0.8).^2 <= 0.01;

%distance
% X_unsafe = y'*y <= 0.3;

X_unsafe = {y'*y <= 0.3, (y-0.7).^2 <= 0.04};


% X_unsafe

%L2 norm
% dist = (x-y)'*(x-y);

%L infinity norm
% dist = [x-y; y-x];

%L1 norm
% dist = {[1; -1]*(x(1)-y(1)), [1; -1]*(x(2)-y(2))}; 

%L3 norm
dist = {[1; -1]*(x(1)-y(1))^3, [1; -1]*(x(2)-y(2))^3}; 


lsupp = unsafe_support(vars);
lsupp.X = X;
% lsupp.X_sys = X_sys;
lsupp.X_init = X_init;
lsupp.param = TH;
lsupp.disturb = W;
lsupp.X_unsafe = X_unsafe;
lsupp.dist = dist;


X0 = lsupp.supp_init();
Xp = lsupp.supp_term();
Xs = lsupp.supp_sys();

f = -x+0.1*(b-0.5);
loc = location_distance(lsupp, f);


%constraint test
% [obj_min, obj_con] = loc.objective_con();
% marg_cons = loc.marg_wass_con(6);
order = 1;
d = 2*order;
[objective, cons_eq, cons_ineq] = loc.all_cons(d)