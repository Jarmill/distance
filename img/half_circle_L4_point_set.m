%Find the point-set L4 distance between a point x and a point y on a
%half-sphere in 3d. This problem is an SOCP


x0 = [0.5; 0.5];
x = sdpvar(2, 1);

y = sdpvar(2, 1);

d = x-y; %distance to minimize

d2 = sdpvar(2, 1); %squared displacement
d4 = sdpvar(2, 1); %fourth power

Cu = [0.25; 0];
Ru = 0.2;

c1f = (-(y(1)-Cu(1))^2 - (y(2)-Cu(2))^2 + Ru^2);
c2f = -y(1)-y(2); 
cons = [c1f; c2f] >= 0;


M12 = cell(2, 1);
M24 = cell(2, 1);
for i = 1:2
    M12{i} = [1, d(i); d(i), d2(i)];
    M24{i} = [1, d2(i); d2(i), d4(i)];
    cons = [cons; M12{i} >= 0; M24{i} >= 0];
end

objective = sum(d4);

opts = sdpsettings('solver', 'mosek', 'savesolverinput', 'true');
vars = [objective; y; d2; d4];
% sol = optimize(cons, objective, opts)
P = optimizer(cons, objective, opts, x, objective)
% sol = P(x0)
P_vec = @(x1, x2) P([x1; x2]);
figure(4)
clf
levellist = linspace(0, 1, 50);
BOX = 1.5;
Fc = fcontour(@(x1,x2) arrayfun(P_vec, x1,x2), [-1, 1, -1, 1]*BOX, 'LevelList', levellist)
