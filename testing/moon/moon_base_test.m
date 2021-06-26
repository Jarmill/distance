h_in = 1.9;
h_out = 2;

%works when h_in, h_out <=1

x_moon = moon_base(h_in, h_out);

figure(700)
clf
hold on
%moon contour
plot(x_moon(1, :), x_moon(2, :), 'k', 'LineWidth', 3)
 


%% test the contour
c_in = [0.5*(1/h_in - h_in)];
r_in = 0.5*(1/h_in + h_in);

c_out = [0.5*(1/h_out - h_out)];
r_out = 0.5*(1/h_out + h_out);

dist = 1.05;

[x_inner, x_outer] = moon_intersection_points(c_in, c_out, dist);


x_cont = moon_contour_base_2(h_in, h_out, dist);

plot(x_cont(1, :), x_cont(2, :), 'r', 'LineWidth', 3)
  scatter(x_inner(1),x_inner(2), 200, 'sk', 'DisplayName', 'Closest Unsafe', 'LineWidth', 2);        
  scatter(x_outer(1),x_outer(2), 200, 'vg', 'DisplayName', 'Closest Unsafe', 'LineWidth', 2);        

   scatter(0,sqrt(dist^2-1), 200, 'ob', 'DisplayName', 'Closest Unsafe', 'LineWidth', 2);        
   scatter(0,-sqrt(dist^2-1), 200, 'ob', 'DisplayName', 'Closest Unsafe', 'LineWidth', 2);        

  
  scatter(0, c_in, 200, 'xk')
  scatter(0, c_out, 200, 'xg')
pbaspect([diff(xlim), diff(ylim), 1])