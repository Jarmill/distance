figure(1)

clf

hold on

% patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none', 'DisplayName', 'Unsafe Set')
   
for i = 1:length(out_sim)
    if i == 1
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
    else
        plot(out_sim{i}.x(:, 1), out_sim{i}.x(:, 2), 'c', 'HandleVisibility', 'Off');
    end
end
plot(out_sim_peak.x(:, 1), out_sim_peak.x(:, 2), 'b', 'DisplayName', 'Closest Traj.', 'LineWidth', 2);       

% scatter(x0_rec(1), x0_rec(2), 200, 'ob', 'DisplayName', 'Closest Initial', 'LineWidth', 2);        
% scatter(xp_rec(1), xp_rec(2), 200, '*b', 'DisplayName', 'Closest Point', 'LineWidth', 2);        
% scatter(xu_rec(1), xu_rec(2), 200, 'sb', 'DisplayName', 'Closest Unsafe', 'LineWidth', 2);        
plot(X0(1, :), X0(2, :), 'k', 'LineWidth', 2, 'DisplayName', 'Initial Set')

% plot([xp_rec(1); xu_rec(1)], [xp_rec(2); xu_rec(2)], ':k', 'DisplayName', 'Closest Distance', 'Linewidth', 1.5)

axis off
 
xlim([-0.2, 1.5])
ylim([-0.6, 0.1])
pbaspect([diff(xlim),diff(ylim),1])