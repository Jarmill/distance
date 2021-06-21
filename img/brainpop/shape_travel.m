%plots of a shape traveling along a trajectory

load('flow_distance.mat', 'out_sim_peak');
shape_color = [179, 0, 255]/255;
side = 0.1;
% R_shape = eye(2);
shape_angle = 5*pi/12;
ang_vel = 1;
R_shape = [cos(shape_angle) -sin(shape_angle); sin(shape_angle) cos(shape_angle)];

 rect_shape_0 = R_shape*side*[1,1,-1,-1,1;-1,1,1,-1,-1];

 
%  Nshape = 5;
 ind_shape = [1; 40; 70; 115];
 figure(1)
 clf
 tl = tiledlayout(1, 2);
%  tl.spacing = 'Compact';
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

FS_title = 16;
 %% translating shape
 ax1 = nexttile;
 hold on
 
   plot(out_sim_peak.x(:, 1), out_sim_peak.x(:, 2), 'b', 'DisplayName', 'Closest Traj.', 'LineWidth', 2);       
   
 for i = 1:length(ind_shape)
     xcurr = out_sim_peak.x(ind_shape(i), :)';
     rect_shape_curr = rect_shape_0 + xcurr;
    patch(rect_shape_curr(1, :), rect_shape_curr(2, :),'k', 'Linewidth', 3,'DisplayName', 'Shape', 'EdgeColor', shape_color, 'FaceColor', 'None')
      scatter(xcurr(1), xcurr(2), 600, '.b', 'DisplayName', 'Closest Shape Center', 'LineWidth', 2);        
        
 end
 
     xlim([-0.4, 1.7])
    ylim([-0.7, 0.3])
    pbaspect([diff(xlim), diff(ylim), 1])
 axis off
 title('Angular Velocity = 0 rad/sec', 'FontSize', FS_title)
 
 
    %% rotating shape
    ax2 = nexttile;
    hold on
    plot(out_sim_peak.x(:, 1), out_sim_peak.x(:, 2), 'b', 'DisplayName', 'Closest Traj.', 'LineWidth', 2);       
   
 for i = 1:length(ind_shape)
     xcurr = out_sim_peak.x(ind_shape(i), :)';
     
     angle_curr = ang_vel * out_sim_peak.t(ind_shape(i));
     
     R_shape_curr = [cos(angle_curr) -sin(angle_curr); sin(angle_curr) cos(angle_curr)];
     
     rect_shape_curr = R_shape_curr*rect_shape_0 + xcurr;
    patch(rect_shape_curr(1, :), rect_shape_curr(2, :),'k', 'Linewidth', 3,'DisplayName', 'Shape', 'EdgeColor', shape_color, 'FaceColor', 'None')
      scatter(xcurr(1), xcurr(2), 600, '.b', 'DisplayName', 'Closest Shape Center', 'LineWidth', 2);        
        
 end
 
     xlim([-0.4, 1.7])
    ylim([-0.7, 0.3])
    pbaspect([diff(xlim), diff(ylim), 1])
    axis off
     title('Angular Velocity = 1 rad/sec', 'FontSize', FS_title)
    
    
 
