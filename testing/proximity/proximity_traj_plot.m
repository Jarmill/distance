T = 10;

% Ax = [0 1; -1 0];
% Ay = [0 -1; 1 -1];

gamma_x = 1;
gamma_y = 2;
m = 2;
van = @(x, gamma) [m*x(2); gamma * (1-(m*x(1)).^2).*(m*x(2))*gamma - m*x(1)];
% van =@(x, gamma) [ 2*x(2,:) ; -0.8*x(1,:) - 10*gamma*(x(1,:).^2-0.21).*x(2,:) ];

fx = @(t, x) van(x, gamma_x);
% fy = @(t, x) van(x, gamma_y);
fy =@(t,x) 2*[x(2); -x(1)];

% fx = @(t, x) Ax * x;
% fy = @(t, y) Ay * y;

x0 = [0; 1];
y0 = [1; 0];
% y0 = [1; 1]/sqrt(2);
% y0 = [0; -1];
% y0 = [0; 0.5];

% [tx_eval, x_eval] = ode15s(fx, [0, T], x0);
% [ty_eval, y_eval] = ode15s(fy, [0, T], y0);
curr_ode_options =   odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep', 0.01);

sol_x = ode15s(fx, [0, T], x0, curr_ode_options );
sol_y = ode15s(fy, [0, T], y0, curr_ode_options );
t_eval = linspace(0, T, 500);

x_eval = deval(sol_x, t_eval);
y_eval = deval(sol_y, t_eval);

c_eval = sum((x_eval - y_eval).^2, 1);

fprintf('Proximity: %0.4f \n', sqrt(min(c_eval)));

figure(1)
clf
hold on
% plot(x_eval(:, 1), x_eval(:, 2));
% plot(y_eval(:, 1), y_eval(:, 2));
% plot3(tx_eval, x_eval(:, 1), x_eval(:, 2));
% plot3(ty_eval, y_eval(:, 1), y_eval(:, 2));
plot3(t_eval, x_eval(1, :), x_eval(2, :));
plot3(t_eval, y_eval(1, :), y_eval(2, :));
view(3)
title('Pair of Trajectories')

% figure(2)
% clf
% % plot(t_eval, c_eval, 'k', 'LineWidth', 3)
% plot(t_eval, sqrt(c_eval), 'k', 'LineWidth', 3)
% title('Proximity between trajectories')
% 
