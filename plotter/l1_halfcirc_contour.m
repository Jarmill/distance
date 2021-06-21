function [x_contour, st_contour]= l1_halfcirc_contour(dist, R, theta_c, Narc)
%     points with an l1 distance of dist to a half-circle with radius R
%     and angle theta_c (the midpoint of the circular arc is at theta_c)

    %handles transformations in case theta is outside the range [0, pi/4]

    %figure out the quadrant and remainder of the tilt angle
%     spins = floor(theta_c / (pi/2));
    
%     theta_c_mod = mod(theta_c, pi/2);
    
    theta_c_mod = wrapTo2Pi(theta_c);
    spins = 0;
    %different than modulus, due to the > sign in the while loop
    while theta_c_mod > pi/2
        theta_c_mod = theta_c_mod - (pi/2);
        spins = spins+1;
    end

    if theta_c_mod > pi/4
        theta_c_base = (pi/2) - theta_c_mod;
        [x_contour, st_contour] = l1_halfcirc_contour_base(dist, R, theta_c_base, Narc);
        x_contour = x_contour([2;1], :);
    else
        theta_c_base = theta_c_mod;
        [x_contour, st_contour] = l1_halfcirc_contour_base(dist, R, theta_c_base, Narc);        
    end
    
    %account for the rotation
    R_spin = [0, -1; 1, 0];
    x_contour = R_spin^(spins) * x_contour;
end
