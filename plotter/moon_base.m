function [outputArg1,outputArg2] = moon_base(inputArg1,inputArg2)
%MOON_BASE Summary of this function goes here
%   Detailed explanation goes here
if h_in == 0
    x_moon_in = [1,-1; 0, 0];
else
    angle_moon_in  = asin(1/r_in);
    if c_in >= 0
        theta_moon_in  = linspace(angle_moon_in, -angle_moon_in, Narc)-pi/2;
    else
        theta_moon_in  = linspace(pi-angle_moon_in, angle_moon_in, Narc)-pi/2;
    end
    x_moon_in = (r_in)*[cos(theta_moon_in); sin(theta_moon_in)] + [0; c_in];
end

angle_moon_out = asin(1/r_out);
theta_moon_out = linspace(-angle_moon_out, angle_moon_out, Narc)-pi/2;

x_moon_out = (r_out)*[cos(theta_moon_out); sin(theta_moon_out)] + [0; c_out];

% x_in_refl = x_in
x_moon = [x_moon_out, x_moon_in];


end

