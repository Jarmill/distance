%try out a tssos decomposition on a problem coming from distance estimation

%distance estimation on the flow system

%decompose the constraint c(x, y) >= w(x) on X times Y for fixed cost c

%This distance constraint does possess term sparsity on "block", but only
%on the first step of the block hierarchy. The second block closure
%operation kills the term sparsity.

%I will need to implement or find the cs-tssos method to hopefully retain
%sparsity while preserving the objective.

x = sdpvar(2,1);
y = sdpvar(2,1);

vars = [x;y];

%L2 distance

c = sum((x-y).^2);

%% indexing powers of a paramterized polynomial of (x, y)
d = 6;
[w, cw] = polynomial(x, d);

%state set 
BOX = 3;

%unsafe set
theta_c = 5*pi/4;       %safe, order 4: 0.2831, order 5: 0.2832
Cu = [0; -0.7];
Ru = 0.5; 

c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 


%lie derivative

%support set in state
All_joint = struct('ineq', [BOX-x.^2; c1f; c2f], 'eq', []);

%% identify monomial powers in support

p_joint = c - w;
data = get_support_monom(p_joint, All_joint, vars, d);

data_old = data;

data = blocks_first(data);

data_first = data;

data = blocks_higher(data);

%this lie derivative is `unblockable', TSSOS with the block setting will
%not decompose this constraint.

%perhaps it could work with other types of problems