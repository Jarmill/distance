%demonstrate that trajectories would collide with the moon if the inner
%circle was flat

% load('moon_result.mat', 'out_sim_multi')
load('moon_result.mat')
%L2 distance bound of  0.1592

figure(54)
clf
hold on

theta = linspace(0, 2*pi, 200);
circ = [cos(theta); sin(theta)];
X0 = C0 + circ*R0;

out_sim = out_sim_multi;
for i = 1:length(out_sim)
    if i == 1
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
    else
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'Off');
    end
end

plot(X0(1, :), X0(2, :), 'k', 'Linewidth', 3, 'DisplayName', 'Initial Set')


%plot the moon
h_in = 0.4;
h_out = 1;
x_moon = moon_base(h_in, h_out);
 

%hugging the curve
moon_center = [0.4;-0.4];
moon_theta = -pi/10;
moon_scale = 0.8;

%moon plot
moon_rot = [cos(moon_theta), sin(-moon_theta); sin(moon_theta), cos(moon_theta)];
x_moon_move = moon_rot*x_moon*moon_scale + moon_center;

x_corner = [-1,1;0,0];
x_corner_move = moon_rot*x_corner*moon_scale + moon_center;


%statistics of the moon
c_in = [0;0.5*(1/h_in - h_in)];
r_in = 0.5*(1/h_in + h_in);

c_out = [0;0.5*(1/h_out - h_out)];
r_out = 0.5*(1/h_out + h_out);


c_in_scale = moon_rot*c_in*moon_scale + moon_center;
c_out_scale = moon_rot*c_out*moon_scale + moon_center;

r_in_scale = moon_scale*r_in;
r_out_scale = moon_scale*r_out;

patch(x_moon_move(1, :), x_moon_move(2, :), 'r', 'EdgeColor', 'None')
plot(x_corner_move(1, :), x_corner_move(2, :), '--r', 'LineWidth', 3)
       xlim([-1, 2.5])
    ylim([-2, 1.5])
    xlabel('x_1')
    ylabel('x_2')
    axis square
