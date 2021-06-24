%(Height, Center-y, Radius) of circles

%Height: minimum y coordinate on circle (as -H)

%inner circle
% h = 0.2;
% h = 0.8;
% h = 1.8;
% h = 2.5;
h = 0.1;
c = 0.5*(1/h - h) ;
r = 0.5*(h + 1/h);

%outer circle
% H = 1;
H = 3;
C = 0.5*(1/H - H);
R = 0.5*(H + 1/H);
 

% dist = 0.5;
% dist = 2.8;
dist = 0.5;
% dist = 1.023
% dist = 2.61;
%when dist > h+c there is a transition, and there will be a pinch point in
%the middle

%% intersection computation
syms x y

eq_in = ((r-dist)^2 == x^2 + (y-c)^2);
eq_side = ((dist)^2 == (x-1)^2 + (y)^2);
eq_out = ((R+dist)^2 == x^2 + (y-C)^2);
[xin, yin] = solve([eq_in,eq_side], [x,y]);
[xout, yout] = solve([eq_out,eq_side], [x,y]);


%% figure out angle ranges
theta_outer = atan2(yout-C, xout);
theta_outer_side = atan2(yout, xout-1);
theta_inner_side = atan2(yin, xin-1);
theta_inner = atan2(yin-c, xin);


Narc = 300;

Nout = floor(Narc * (theta_outer+pi/2)/(2*pi))+1;

Nside = floor(Narc * (theta_outer_side - theta_inner_side)/(2*pi))+1;

Nin = floor(Narc * (theta_inner+pi/2)/(2*pi))+1;

theta_out = linspace(-pi/2, theta_outer, Nout);
theta_side = linspace(theta_outer_side,wrapTo2Pi(theta_inner_side), Nout);
theta_in = linspace(theta_inner, -pi/2, Nout);


x_out = (R+dist)*[cos(theta_out); sin(theta_out)] + [0; C];
x_side = (dist)*[cos(theta_side); sin(theta_side)] + [1; 0];
x_in = (r-dist)*[cos(theta_in); sin(theta_in)] + [0; c];

x_half = [x_out, x_side, x_in];

x_refl = diag([-1,1])*x_half(:, end:-1:1);

x_moon = [x_half, x_refl];

% angle_low = 


%% plot for visualization
Nth = 201;
theta = linspace(0, 2*pi, Nth); 
circ = [cos(theta); sin(theta)];

%circles forming the moon
circ_inner = r*circ + [0; c];
circ_outer = R*circ + [0; C];

%circles forming the distance contours
% circ_inner_dist = sqrt(r^2 - dist)*circ + [0; c];
% circ_outer_dist = sqrt(R^2 + dist)*circ + [0; C];
% circ_left_dist  = sqrt(dist)*circ + [-1; 0];
% circ_right_dist = sqrt(dist)*circ + [1; 0];

circ_inner_dist = (r - dist)*circ + [0; c];
circ_outer_dist = (R + dist)*circ + [0; C];
circ_left_dist  = (dist)*circ + [-1; 0];
circ_right_dist = (dist)*circ + [1; 0];


figure(35)
clf
hold on
plot(circ_inner(1, :), circ_inner(2, :), 'k')
plot(circ_outer(1, :), circ_outer(2, :), 'k')
scatter(0, C, 600, '.g')
scatter(0, c, 600, '.r')

%circles
plot(circ_inner_dist(1, :), circ_inner_dist(2, :), 'r')
plot(circ_outer_dist(1, :), circ_outer_dist(2, :), 'g')
plot(circ_left_dist(1, :), circ_left_dist(2, :), 'b')
plot(circ_right_dist(1, :), circ_right_dist(2, :), 'b')

%moon contour
plot(x_moon(1, :), x_moon(2, :), 'm', 'LineWidth', 3)

scatter([1,-1], [0,0], 200, 'sk')
scatter([xin,-xin], [yin,yin], 200, 'hk')
scatter([xout,-xout], [yout,yout], 200, 'vk')

pbaspect([diff(xlim), diff(ylim), 1])

