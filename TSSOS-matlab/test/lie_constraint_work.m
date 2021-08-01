t = sdpvar(1,1);
x = sdpvar(1,1);
y = sdpvar(1,1);

vars = [t; x;y];

%% indexing powers of a paramterized polynomial of (x, y)
n = 2;
d = 8;
[v, cv] = polynomial(vars, d);

%dynamics of van der pol oscillator
f = [1; 2*y; -0.5*x - 10*(x^2-0.2)*y];

%lie derivative
Lv = jacobian(v, vars)*f;

%support set in state
X = struct('ineq', [t*(1-t); 1-x^2; 1-y^2], 'eq', []);

%% identify monomial powers in support
%this line finds all monomial powers in Lv
% [cc_Lv, mm_Lv] = coefficients(Lv, [x; y]);
[ex_Lv, base_Lv] = getexponentbase(Lv, vars);
ex_Lv_unique = unique(full(ex_Lv), 'rows', 'sorted');
%[exponents,base]=getexponentbase(p,x)


data = get_support_monom(Lv, X, vars, d);

data_old = data;

data = blocks_first(data);


%this lie derivative is `unblockable', TSSOS with the block setting will
%not decompose this constraint.

%perhaps it could work with other types of problems