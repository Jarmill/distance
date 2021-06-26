% load('moon_result.mat', 'out_sim_multi')
load('moon_result.mat')
%L2 distance bound of  0.1592

figure(53)
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

%same coordinates as half-circle example
% moon_center = [0;-0.7];
% moon_theta = -pi/4;
% moon_scale = 0.5

%moon plot
moon_rot = [cos(moon_theta), sin(-moon_theta); sin(moon_theta), cos(moon_theta)];
x_moon_move = moon_rot*x_moon*moon_scale + moon_center;

%distance contour plot
xdist_raw = moon_contour_base_2(h_in, h_out, dist_rec/moon_scale);

xdist_move = moon_rot*xdist_raw*moon_scale + moon_center;


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
 

patch(xdist_move(1, :), xdist_move(2, :),'r', 'Linewidth', 3,'DisplayName', 'L1 Contour', ...        
    'EdgeColor', 'r', 'FaceColor', 'None')

f = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];

[t_opt, x_opt] = ode45(f, [0, 5], mom_rec.x0);
out_sim_peak.t = t_opt;
out_sim_peak.x = x_opt;

    plot(out_sim_peak.x(:, 1), out_sim_peak.x(:, 2), 'b', 'DisplayName', 'Closest Traj.', 'LineWidth', 2);       

    scatter(mom_rec.x0(1), mom_rec.x0(2), 200, 'ob', 'DisplayName', 'Closest Initial', 'LineWidth', 2);        
    scatter(mom_rec.xp(1), mom_rec.xp(2), 200, '*b', 'DisplayName', 'Closest Point', 'LineWidth', 2);        
    scatter(mom_rec.y(1), mom_rec.y(2), 200, 'sb', 'DisplayName', 'Closest Unsafe', 'LineWidth', 2);        

    plot([mom_rec.xp(1); mom_rec.y(1)], [mom_rec.xp(2); mom_rec.y(2)], ':k', 'DisplayName', 'Closest Distance', 'Linewidth', 1.5)

       xlim([-1, 2.5])
    ylim([-2, 1.5])
    xlabel('x_1')
    ylabel('x_2')
    axis square
