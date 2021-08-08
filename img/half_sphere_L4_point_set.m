%Find the point-set L4 distance between a point x and a point y on a
%half-sphere in 3d. This problem is an SOCP
x = sdpvar(3, 1);

y = sdpvar(3, 1);

d = x-y; %distance to minimize

d2 = sdpvar(3, 1); %squared displacement
d4 = sdpvar(3, 1); %fourth power

Cu = [0.25; 0; 0];
Ru = 0.2;

c1f = (-(y(1)-Cu(1))^2 - (y(2)-Cu(2))^2 - (y(3)-Cu(3))^2 + Ru^2);
c2f = -y(3); 
cons = [c1f; c2f] >= 0;


M12 = cell(3, 1);
M24 = cell(3, 1);
for i = 1:3
    M12{i} = [1, d(i); d(i), d2(i)];
    M24{i} = [1, d2(i); d2(i), d4(i)];
    cons = [cons; M12{i} >= 0; M24{i} >= 0];
end

objective = sum(d4);

opts = sdpsettings('solver', 'mosek');
vars = [objective; y; d2; d4];
% sol = optimize(cons, objective, opts)
P = optimizer(cons, objective, opts, x, objective)
% sol = P(x0)
P_vec = @(x1, x2, x3) P([x1; x2; x3]);
figure(4)
clf
% levellist = linspace(0, 1, 50);
% levellist =  0.041115989882192;

%L4^4 bound
levellist =  2.857873504325048e-06;

P_vec_bound = @(x1, x2, x3) P([x1; x2; x3]) - levellist;

% BOX = 1.5;
figure(1)
falpha = 0.4;
Fc = fimplicit3(@(x1,x2, x3) arrayfun(P_vec_bound, x1,x2, x3), [-0.25, 0.75, -0.5, 0.5, -0.5, 0.25], 'EdgeColor', 'None', 'FaceColor', 'r', 'FaceAlpha', falpha)
% Fc = fcontour3(@(x1,x2, x3) arrayfun(P_vec, x1,x2, x3), [-0.25, 0.75, -0.5, 0.5, -0.5, 0.25], 'LevelList', levellist)
