%Find the point-set L4 distance between a point x and a point y on a
%half-sphere in 3d. This problem is an SOCP


x0 = [0.5; 0.5; 0.5];
x = sdpvar(3, 1);

y = sdpvar(3, 1);

d = x-y; %distance to minimize

d2 = sdpvar(3, 1); %squared displacement
d4 = sdpvar(3, 1); %fourth power

Cu = [0.25; 0; 0];
Ru = 0.2;

c1f = (-(y(1)-Cu(1))^2 - (y(2)-Cu(2))^2  - (y(3) - Cu(3))^2 + Ru^2);
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

opts = sdpsettings('solver', 'mosek', 'savesolverinput', 'true');
vars = [objective; y; d2; d4];
% sol = optimize(cons, objective, opts)
P = optimizer(cons, objective, opts, x, vars)
sol = P(x0)