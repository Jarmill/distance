%contour plot of l1 distance of a circle
Nth = 10000;
theta = linspace(0, pi/2, Nth);
R = 1;

% https://math.stackexchange.com/questions/2202124/l1-distance-to-l2-sphere

circ = [cos(theta); sin(theta)];

x = [3; 1];

[dist, x_circ] = l1_circ_emp(x, R, circ);

PLOT = 0;

emp_handle = @(x, y) l1_circ_emp([x; y], R, circ);

quart_handle = @(x, y) l1_circ_quarter([x; y], R);

if PLOT
% figure(1)
% clf
% level_list = linspace(0, 4, 30);

% 
% fsurf(l1_handle, [0, 3, 0, 3])
% figure(2)
% 
% clf
% hold on 
% fcontour(l1_handle, [0, 3, 0, 3],...
% 'LevelList', level_list)
% plot(R*circ(1, :), R*circ(2, :), 'k')
% axis square


Ngrid = 100;
[XX, YY] = meshgrid(linspace(0, 3, Ngrid));

l1_quart = arrayfun(quart_handle, XX, YY);
emp_quart = arrayfun(emp_handle, XX, YY);

% surf(XX, YY, l1_quart-emp_quart)
contour(XX, YY, l1_quart)
end

x_test = [3; 0.45];
[d_emp, x_emp] = emp_handle(x_test(1), x_test(2));
[d_quart, x_quart] = quart_handle(x_test(1), x_test(2));

figure(3)
clf
hold on
Ndist = 21;
dist_range = linspace(0, 4, Ndist);
for i = 1:Ndist
    x_contour = l1_circ_quart_contour(dist_range(i), R, 100);
    plot(x_contour(1, :), x_contour(2, :), 'k')
    axis square
end
function [dist, x_circ]= l1_circ_emp(x, R, circ)
    %empirical distance based on samples 
    if norm(x) < R
        dist = 0;
        x_circ = x;
    else    
        l1_emp = sum(abs(x-R*circ), 1);

        [dist, ind_opt] = min(l1_emp);

        x_circ = R*circ(:, ind_opt);
    end
end

function [dist, x_circ] = l1_circ_quarter(x, R)
    %distance from point to quarter-circle
    x_crit = R/sqrt(2)*[1; 1];
    if norm(x) < R
        dist = 0;
        x_circ = x;
    else
        if all(x > x_crit)
            %top right quarter-space starting from x_crit            
            x_circ = x_crit;
        else
            if x(1) > x_crit
                %horizontal dominance
                x_circ = [sqrt(R^2 - x(2)^2); x(2)];
            else
                %vertical dominance
                x_circ = [x(1); sqrt(R^2 - x(1)^2)]; 
            end
        end    
        dist = sum(abs(x - x_circ));
    end

end

function x_contour = l1_circ_quart_contour(dist, R, Narc)
%     points with an l1 distance of dist to a circle with radius R
%     Narcsamples of circle arcs
%     theta_c = atan2(R/sqrt(2), R + dist);
%     theta_c = atan2(R/sqrt(2), R/sqrt(2) + dist);
    theta_c = pi/4;
%     theta_c = 2*pi;

    theta = linspace(0, theta_c, Narc);
    arc = R*[cos(theta); sin(theta)] + [dist; 0];
    arc_refl = arc([2, 1], end:-1:1);
    %reflect over the arc
%     u = [-1; 1]/sqrt(2);
%     H = eye(2) - 2*u*u'
    x_contour = [arc, arc_refl];    
end