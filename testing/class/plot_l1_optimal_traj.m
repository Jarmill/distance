load('l1_opt_result.mat')
%L1 distance bound of 0.4003

% out_sim = opt_result.out_sim;
% out_sim_peak = opt_result.out_sim_peak;

% smp = struct('x', @() sphere_sample(1, 2)'*R0 + C0);

% LS = sampler_base(PM.loc, smp);

% Ntraj = 50;
% out_sim_multi = LS.sample_traj_multi(Ntraj, lsupp.Tmax);

% out_sim_peak = LS.sample_traj(0, s_init.x, 5); 


figure(53)
clf
hold on

theta = linspace(0, 2*pi, 200);
circ = [cos(theta); sin(theta)];
X0 = C0 + circ*R0;

theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
circ_half = [cos(theta_half_range); sin(theta_half_range)];
Xu = Cu + circ_half* Ru;

out_sim = out_sim_multi;
for i = 1:length(out_sim)
    if i == 1
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
    else
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'Off');
    end
end

plot(X0(1, :), X0(2, :), 'k', 'Linewidth', 3, 'DisplayName', 'Initial Set')
patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')

dist = sol.obj_rec;
x_dist = l1_halfcirc_contour(sol.obj_rec, Ru, theta_c, 500) + Cu;

patch(x_dist(1, :), x_dist(2, :),'r', 'Linewidth', 3,'DisplayName', 'L1 Contour', 'EdgeColor', 'r', 'FaceColor', 'None')       

% if optimal_pt
    plot(out_sim_peak.x(:, 1), out_sim_peak.x(:, 2), 'b', 'DisplayName', 'Closest Traj.', 'LineWidth', 2);       

    scatter(opt_result.x0(1), opt_result.x0(2), 200, 'ob', 'DisplayName', 'Closest Initial', 'LineWidth', 2);        
    scatter(opt_result.xp(1), opt_result.xp(2), 200, '*b', 'DisplayName', 'Closest Point', 'LineWidth', 2);        
    scatter(opt_result.y(1), opt_result.y(2), 200, 'sb', 'DisplayName', 'Closest Unsafe', 'LineWidth', 2);        

    plot([opt_result.xp(1); opt_result.y(1)], [opt_result.xp(2); opt_result.y(2)], ':k', 'DisplayName', 'Closest Distance', 'Linewidth', 1.5)
% end
       xlim([-1, 2.5])
    ylim([-2, 1.5])
    xlabel('x_1')
    ylabel('x_2')
    axis square
    
%         if status_feas == 0
        %feasible program 
%         title_str = (['Unsafe, false L_2 distance bound is ', num2str(dist_rec, 3)]);
%     else    
%         title_str = (['L_1 distance bound is ', num2str(dist_rec, 3)]);        
%     end
%     title(title_str, 'FontSize' , FS_title)
   
