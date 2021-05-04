load('occ_data.mat')

%should have put this in the structure
T = 2;
xrange = [-4, 4];


ind_special = 2000;
probestride = 20; %=40
x_probe = [1, probestride:probestride:size(x_traj, 1)];


%% compute the boxes
% Red: neither blue or any cyan passes through
% Green: some cyan passes through, not blue
% Black: some cyan and blue pass through
t_red = [1, 1.6];
x_red = [2, 3];

t_green = [0.25, 0.8];
x_green =  [-1, 2];

t_black = [0.5, 1.25];
x_black = [-3, -1.5];


%% trajectories
figure(1)
clf
hold on
plot(t_traj, x_traj(x_probe, :), 'c')
plot(t_traj, x_traj(ind_special, :), 'b', 'LineWidth', 3)

%% boxes
patch(t_red([1,1,2,2,1]), x_red([2,1,1,2,2]), 'r', ...
    'EdgeColor', 'g', 'LineWidth', 3,  'FaceColor', 'None')

patch(t_green([1,1,2,2,1]), x_green([2,1,1,2,2]), 'r', ...
    'EdgeColor', 'r', 'LineWidth', 3,  'FaceColor', 'None')


patch(t_black([1,1,2,2,1]), x_black([2,1,1,2,2]), 'r', ...
    'EdgeColor', 'k', 'LineWidth', 3,  'FaceColor', 'None')

xlabel('t')
ylabel('x(t)')
% title('Distributions', 'FontSize', 16)