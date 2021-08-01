%% test with indexing a single polynomial of (x, y)
x = sdpvar(1,1);
y = sdpvar(1,1);

m = x^2 * y^3;

p = m + 2*x*y - 3*y + 1;

% [cc, mm] = coefficients(p);
[ex, base] = getexponentbase(p, [x; y]);

[ex_y, base_y] = getexponentbase(p, [y]);

%cc:    coefficients of polynomial
%mm:    monomials of polynomial associated with coefficients

%need to find a way to extract out exponents of monomials

