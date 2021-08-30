%estimate distance between two trajectoires of a linear system using an
%LMI. Then scale this up to sets of initial conditions and more complicated
%systems.

%% information about dynamics
SOLVE = 0;
PLOT = 1;

Tmax = 10;

Ax = [0 1; -1 0];
Ay = [0 -1; 1 -1];

fx = @(t, x) Ax * x;
fy = @(t, y) Ay * y;

X0 = [0; 1];
Y0 = [1; 0];

BOX = 1;

order = 2;
d = 2*order;
if SOLVE
mset clear
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'))
%% set up measures
mpol('tx0', 1, 1);
mpol('x0', 2, 1);
mu0x = meas([tx0; x0]);

mpol('ty0', 1, 1);
mpol('y0', 2, 1);
mu0y = meas([tx0; x0]);

mpol('tx_occ', 1, 1);
mpol('x_occ', 2, 1);
mu_occx = meas([tx_occ; x_occ]);

mpol('ty_occ', 1, 1);
mpol('y_occ', 2, 1);
mu_occy = meas([ty_occ; y_occ]);


mpol('txp', 1, 1);
mpol('xp', 2, 1);
mupx = meas([txp; xp]);

mpol('typ', 1, 1);
mpol('yp', 2, 1);
mupy = meas([typ; yp]);


%wasserstein
mpol('xw', 2, 1)
mpol('yw', 2, 1)
eta = meas([xw; yw]);

%% support constraints
X0_supp = x0 == X0;
Y0_supp = y0 == Y0;


supp_con = [tx0 == 0; ty0 == 0; tx_occ*(1-tx_occ)>=0; ty_occ*(1-ty_occ)>=0; 
    txp*(1-txp) >=0; typ*(1-typ) >=0;
    x_occ.^2 <= BOX^2; xp.^2 <= BOX^2; 
    y_occ.^2 <= BOX^2; yp.^2 <= BOX^2; 
    X0_supp; Y0_supp;
    xw.^2 <= BOX^2; yw.^2 <= BOX^2;    
    ];



%% moment constraints

%liouville constraint
mom_init_x = mom(mmon([tx0; x0], d));
mom_peak_x = mom(mmon([txp; xp], d));

mom_init_y = mom(mmon([ty0; y0], d));
mom_peak_y = mom(mmon([typ; yp], d));


MOM_SUBS =1;

vx  = mmon([tx_occ; x_occ], d);
vy  = mmon([ty_occ; y_occ], d);
% f_occ = subs_vars(f, x, x_occ);
fx_occ = fx(tx_occ, x_occ);
fy_occ = fy(ty_occ, y_occ);

Ax_occ = mom(diff(vx, tx_occ) + diff(vx, x_occ)*fx_occ*Tmax); 
Ay_occ = mom(diff(vy, ty_occ) + diff(vy, y_occ)*fy_occ*Tmax); 

if MOM_SUBS
    Liou_con_x =( mom_peak_x == mom_init_x + Ax_occ);
    Liou_con_y =( mom_peak_y == mom_init_y + Ay_occ);
    
    Liou_con = [Liou_con_x; Liou_con_y];
else
    Liou_x = Ax_occ + (mom_init_x - mom_peak_x);
    Liou_y = Ay_occ + (mom_init_y - mom_peak_y);
    Liou_con = ([Liou_x; Liou_y] == 0);
end

%marginal between peak and wass

% mom_wass_x = mom(mmon(xw), d);
% mom_wass_y = mom(mmon(yw), d);

if MOM_SUBS
    Time_con = mom(mmon(txp, d))==mom(mmon(typ,d));
else
    Time_con = (mom(mmon(txp,d))-mom(mmon(typ,d))) == 0;
end


% ypx = mom(mmon(xp, d));
% ywx = mom(mmon(xw, d));
if MOM_SUBS
%     Wmarg_con = (ywx==ypx); 

    Wmarg_con = [mom(mmon(xp, d))==mom(mmon(xw, d)); mom(mmon(yp, d))==mom(mmon(yw, d))];
%     mass_con = mass(mu0)==1;
else
    Wmarg = [mom(mmon(xp, d))-mom(mmon(xw, d)); mom(mmon(yp, d))-mom(mmon(yw, d))];
    Wmarg_con = Wmarg==0;   
%     mass_con = mass(mu0)-1=0;
end
mom_con = [Liou_con; Wmarg_con; Time_con; mass(mu0x)==1];


%objective
mom_dist = mom(sum((xw-yw).^2));
objective = min(mom_dist);

%% solve problem
%Input LMI moment problem
P = msdp(objective, ...
    mom_con, supp_con);

%solve LMIP moment problem
[status, obj, m, dual_rec] = msol(P);    

%% analyze solutions

dist_rec = sqrt(double(mom_dist));
disp(['distance bound: ', num2str(dist_rec)])

%capture the moment matrices
M0x = double(mmat(mu0x));
M0y = double(mmat(mu0y));
Mpx = double(mmat(mupx));
Mpy = double(mmat(mupy));
Meta = double(mmat(eta));
Moccx = double(mmat(mu_occx));
Moccy = double(mmat(mu_occy));

%save time by only computing moments of corners if necessary
n = 2;
M0x_1 = M0x(1:(n+2), 1:(n+2));
Mpx_1 = Mpx(1:(n+2), 1:(n+2));
M0y_1 = M0y(1:(n+2), 1:(n+2));
Mpy_1 = Mpy(1:(n+2), 1:(n+2));
Meta_1 = Meta(1:(2*n+1), 1:(2*n+1));
rankpx = rank(Mpx_1, 1e-3);
rank0x = rank(M0x_1, 1e-3);
rankpy = rank(Mpy_1, 1e-3);
rank0y = rank(M0y_1, 1e-3);
ranketa = rank(Meta_1, 1e-3);

%extract putative optima if measures have low-rank moment matrices
xp_rec = double(mom(xp));
x0_rec = double(mom(x0));
yp_rec = double(mom(yp));
y0_rec = double(mom(y0));
tp_rec = Tmax*double(mom(txp));


optimal_pt = all([rankpx; rank0x; rankpy; rank0y; ranketa]==1);
end

if PLOT
    curr_ode_options =   odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'MaxStep', 0.01);

    sol_x = ode15s(fx, [0, T], x0, curr_ode_options );
    sol_y = ode15s(fy, [0, T], y0, curr_ode_options );
    t_eval = linspace(0, T, 100);

    x_eval = deval(sol_x, t_eval);
    y_eval = deval(sol_y, t_eval);

    c_eval = sum((x_eval - y_eval).^2, 1);

    
    figure(1)
clf
hold on
% plot(x_eval(:, 1), x_eval(:, 2));
% plot(y_eval(:, 1), y_eval(:, 2));
% plot3(tx_eval, x_eval(:, 1), x_eval(:, 2));
% plot3(ty_eval, y_eval(:, 1), y_eval(:, 2));
plot3(t_eval, x_eval(1, :), x_eval(2, :));
plot3(t_eval, y_eval(1, :), y_eval(2, :)); 
if optimal_pt
    scatter3([0, 0], [x0_rec(1); y0_rec(1)], [x0_rec(2); y0_rec(2)], 200, 'ok')
    plot3([1, 1]*tp_rec, [xp_rec(1); yp_rec(1)], [xp_rec(2); yp_rec(2)], '--*k', 'MarkerSize', 15)
%     plot3(
end
view(3)
title('Pair of Trajectories')

figure(2)
clf
% plot(t_eval, c_eval, 'k', 'LineWidth', 3)
hold on
plot(t_eval, sqrt(c_eval), 'k', 'LineWidth', 3)
plot([0, Tmax], dist_rec*[1, 1], '--r', 'LineWidth', 3);
title('Proximity between trajectories')


end