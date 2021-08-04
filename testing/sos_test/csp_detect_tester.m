%try out a tssos decomposition on a problem coming from distance estimation

%distance estimation on the flow system

%decompose the constraint c(x, y) >= w(x) on X times Y for fixed cost c

%This distance constraint does possess term sparsity on "block", but only
%on the first step of the block hierarchy. The second block closure
%operation kills the term sparsity.

%I will need to implement or find the cs-tssos method to hopefully retain
%sparsity while preserving the objective.

x = sdpvar(4,1);
y = sdpvar(4,1);

vars = [x;y];

%L2 distance

c = sum((x-y).^2);

%% indexing powers of a paramterized polynomial of (x, y)
d = 4;
[w, cw] = polynomial(x, d);

%state set 
BOX = 3;

%unsafe set
theta_c = 5*pi/4;       %safe, order 4: 0.2831, order 5: 0.2832
Cu = [0; -0.7];
Ru = 0.5; 



%lie derivative

%support set in state
All_joint = struct('ineq', [BOX-x.^2; 1-sum(y.^2); y(1)], 'eq', []);

%% identify monomial powers in support

p_joint = c - w;

[cliques] = find_csp(p_joint, All_joint, vars)

% data_first = data;
% 
% data = blocks_higher(data);

%this lie derivative is `unblockable', TSSOS with the block setting will
%not decompose this constraint.

%perhaps it could work with other types of problems