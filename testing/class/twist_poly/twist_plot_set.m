%twist dynamics
A_true = [-1 1 1; -1 0 -1; 0 1 -2];
B_true = [-1 0 -1;
          0 1 1;
          1 1 0]/2;

f_true = @(t,x) A_true*x - B_true*(4*x.^3 - 3*x);

%unsafe set

Cu = [0.25; 0; 0];
Ru = 0.2;
[Xs, Ys, Zs] = sphere(40);
Xu = Ru*Xs + Cu(1);
Yu = Ru*Ys + Cu(2);
Zu = Ru*Zs + Cu(3);

% sampling
% [tt, xx] 
C0 = [-0.5; 0; 0];
R0 = 0.2;
% sampler_x = @() [ball_sample(1,2), 0]'*R0 + C0;
sampler_x = @() [sphere_sample(1,3)]'*R0 + C0;


% Tmax = 1000;
Tmax = 5;

Nsample = 150;
out_sim = cell(Nsample, 1);
figure(1)
clf
hold on
for i = 1:Nsample
    x0 = sampler_x();
% [tt, xt] = ode45(f_true, [0, Tmax], C0);
[tt, xt] = ode45(f_true, [0, Tmax], x0);
out_sim{i} = struct;
out_sim{i}.t = tt;
out_sim{i}.x = xt;
% [min(xt, [], 1);
% max(xt, [], 1)]
plot3(xt(:, 1), xt(:, 2), xt(:, 3))
scatter3(xt(1, 1), xt(1, 2), xt(1, 3), 100, 'k')
end

FS_axis = 12;
xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', FS_axis);
ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', FS_axis);
zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', FS_axis);


surf(Xu, Yu, Zu)
view(3)