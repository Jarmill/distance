
rng(30, 'twister')
% lam = [-0.5 5; -5 -0.5];
% ba = randn(2,2);

% A = ba'*lam*ba;
A = [-0.3, -4.3; 3.7, -0.8]

x0 = [2; 0];
T = 5;
tspan = linspace(0, T, 400);
[t, x] = ode45(@(t, x) A*x, tspan, x0);

xu = [0; 0.5];
dist = sqrt(sum((x'-xu).^2, 1)');
[mindist, idist] = min(dist);
xdist = x(idist, :);

%perform the plot
figure(1)
clf
tiledlayout(3,1)
nexttile([2, 1])
hold on
plot(x(:, 1), x(:, 2))
scatter(xu(1), xu(2), 400, '.k')
scatter(x0(1), x0(2), 200, 'ok')
plot([xu(1); xdist(1)], [xu(2), xdist(2)], '--b')
scatter(xdist(1), xdist(2), 200, '*b')
title('Trajectory Plot')
xlabel('x_1')
ylabel('x_2')

nexttile
hold on
plot(tspan, dist)
scatter(tspan(idist), dist(idist), 200, '*b')
title('Distance to Unsafe Set')
xlabel('time')
ylabel('L2 distance')