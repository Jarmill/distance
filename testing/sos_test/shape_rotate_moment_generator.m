%generate moment and localizing matrix for rotating shape system


syms x y c s sx sy;
% vars = [c; s; sx; sy; x; y];
vars = [c, s, sx, sy, x, y];
%(x, y) center of shape
%(c, s) cosine and sine of rotating angle 
%(sx, sy) local coordinates on shape

BOX = 3;
shape_size = 0.5;

%inequality constraints of state
g = [BOX^2 - [x; y].^2; shape_size^2 - [sx; sy].^2];

%equality constraint, grobner basis
h = c^2 + s^2 -1; %==0

%transformation
shape_transform = [x; y] + [c, -s; s, c]*[sx; sy];

%% compute polynomials

%pushforward operation in the rotation doubles the system order
order = 4;
order2 = 2*order;
d = 2*order2;

Nbin = 1;
Nstd = 5;

%fbasis are monomials on rows and columns of the moment matrix
[fbasis, index] = momentPowers(Nbin, Nstd, order2);

%monomials on rows
monom_row = prod((ones(length(index), 1)*vars).^fbasis, 2);

m = size(fbasis, 1)
temp = bsxfun(@plus,kron(ones(m,1),fbasis),kron(fbasis,ones(m,1)));

%fmom are monomials inside the moment matrix
[fmom, ia, ic] = unique(temp, 'rows');

%elements that require polynomial reduction
ind_reduce = (fmom(:, 1) == 2);
% n_reduce = sum(ind_reduce);
% n_std = length(ia) - n_reduce;


monom_inner = prod((ones(length(ia), 1)*vars).^fmom, 2);

monom_reduce_pre = monom_inner(ind_reduce);

% monom_reduce = poly_reduce_vec(monom_reduce_pre, h, vars);

%this reduction rule is very easy
monom_reduce = (monom_reduce_pre/(c^2))*(1 - s^2);


monom_inner(ind_reduce) = monom_reduce;

M = reshape(monom_inner(ic), m, m);


