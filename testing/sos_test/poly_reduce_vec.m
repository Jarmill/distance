function [r_vec] = poly_reduce_vec(p_vec,g, vars)
%POLY_REDUCE_VEC reduce a vector of polynomials (p_vec) by the ideal formed
%by g
%   Detailed explanation goes here
    if nargin == 2
        r_vec = arrayfun(@(m) polynomialReduce(m, g), p_vec);
    else
        r_vec = arrayfun(@(m) polynomialReduce(m, g, vars), p_vec);
    end
end

