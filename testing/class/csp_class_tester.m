mset clear
n = 4;
mpol('rr', 3, 1)
mpol('x', n, 1);
mpol('y', n, 1);

BOX = 3;
R = 1;
supp_X = BOX^2-x.^2 >= 0;
supp_Xu = [y(1) >= 0; sum(y.^2)<=R^2];
vars = struct('x', x, 'y', y);
mcsp = meas_wass_csp(vars, supp_X, supp_Xu);

d = 3;
mx = mcsp.mom_monom_x(d);
my = mcsp.mom_monom_y(d);

con = mcsp.overlap_con(d);

c = sum((x-y).^2);
c1 = sum((x(1)-y(1)).^2);

% eval(c, [x; y], zeros(2*n, 1))
% eval(c1, [x; y], zeros(2*n, 1))

% p_public = get_representation(c1, [x; y]);

mom_c = mcsp.mom_objective(c);