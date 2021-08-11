%Find the set-set L2 distance between a point x on a rotating square 
%and a point y on a rotated half-circle. Evaluate this distahce along
%trajectories


shape_size = 0.1;

shape_angle = 5*pi/12;
R_shape = [cos(shape_angle) -sin(shape_angle); sin(shape_angle) cos(shape_angle)];
shape_size = 0.1;


x0 = [0.5; 0.5];
% x = x0;

x = sdpvar(2, 1);
% x = [1; 1]
s = sdpvar(2, 1);
y = sdpvar(2, 1);


shape_transform = x + R_shape*s*shape_size;

d = shape_transform - y;
% d = x+s-y; %distance to minimize

% d2 = sdpvar(2, 1); %squared displacement
% d4 = sdpvar(2, 1); %fourth power

Cu = [0; -0.7];
%Cu = [2.5; 0];
Ru = 0.5; 
theta_c = 5*pi/4;  
c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 

cons = [c1f; c2f] >= 0;
cons = [cons; s >= -1; s <= 1];

% M12 = cell(3, 1);
% M24 = cell(3, 1);
% for i = 1:3
%     M12{i} = [1, d(i); d(i), d2(i)];
%     M24{i} = [1, d2(i); d2(i), d4(i)];
%     cons = [cons; M12{i} >= 0; M24{i} >= 0];
% end

% objective = sum(d4);
objective = norm(d, 2);

opts = sdpsettings('solver', 'mosek', 'savesolverinput', 'true');
% vars = [objective; y; d2; d4];
% sol = optimize(cons, objective, opts)
vars = [objective; y; s];
P = optimizer(cons, objective, opts, x, vars);

%% evaluate set-set distances
load('flow_distance.mat', 'out_sim_peak');
out_sim_peak.dist_shape = zeros(length(out_sim_peak.t),1);
for i = 1:length(out_sim_peak.t)
     xcurr = out_sim_peak.x(i, :)';
     Pcurr = P(xcurr);
     out_sim_peak.dist_shape(i) = Pcurr(1);
%      angle_curr = ang_vel * out_sim_peak.t(ind_shape(i));
     
%      R_shape_curr = [cos(angle_curr) -sin(angle_curr); sin(angle_curr) cos(angle_curr)];
     
end

% sol = P(x0)