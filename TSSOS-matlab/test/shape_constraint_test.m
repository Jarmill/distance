%try out a tssos decomposition on a problem coming from distance estimation

%distance estimation with shapes

%decompose the constraint w(A(s; omega)) >= z(omega) on S times Omega

%This distance constraint does possess term sparsity on "block", but only
%on the first step of the block hierarchy. The second block closure
%operation kills the term sparsity.

%I will need to implement or find the cs-tssos method to hopefully retain
%sparsity while preserving the objective.

omega = sdpvar(2,1);
s = sdpvar(2,1);
x = sdpvar(2,1);

vars = [omega;s];

%L2 distance

%coordinte transformation
% A = s + omega;

%% indexing powers of a paramterized polynomial of (x, y)
d = 6;
[w, cw] = polynomial(x, d);
[z, cz] = polynomial(omega, d);

% %state set 
% BOX = 3;
% 
% %unsafe set
% theta_c = 5*pi/4;       %safe, order 4: 0.2831, order 5: 0.2832
% Cu = [0; -0.7];
% Ru = 0.5; 
% 
% c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;
% w_c = [cos(theta_c); sin(theta_c)];
% c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 
% 
% 
% %lie derivative
% 
% %support set in state
% All_joint = struct('ineq', [BOX-x.^2; c1f; c2f], 'eq', []);

BOX = 3;

shape_angle = 5*pi/12;
shape_side = 0.1;

S_Omega = struct('ineq', [BOX^2-omega(1:2).^2; shape_side^2 - s.^2], 'eq', []);


%coordinate transformation
R_shape = [cos(shape_angle) -sin(shape_angle); sin(shape_angle) cos(shape_angle)];
shape_transform = omega(1:2) + R_shape*s;

%% identify monomial powers in support
w_A = replace(w, x, shape_transform);
p_shape = w_A - z;
data = get_support_monom(p_shape, S_Omega, vars, d);

data_old = data;

data = blocks_first(data);
% 
data_first = data;
% 
% data = blocks_higher(data);

%this lie derivative is `unblockable', TSSOS with the block setting will
%not decompose this constraint.

%perhaps it could work with other types of problems