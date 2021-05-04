%visualizes slices of an occupation measure
%inspired by Figure 1 of https://groups.csail.mit.edu/robotics-center/public_papers/Majumdar14.pdf

%Author: Jared Miller, May 3 2021

SAMPLE = 0;
KDE = 1;
PLOT_TRAJ = 0;
PLOT_OCC = 0;

f = @(t, x) -x^3 - x^2 + 2*x;

T = 2;

if SAMPLE
    xrange = [-4, 4];
    
    Xsample = 20000;
    Tsample = 200;
       
    options = odeset('AbsTol', 1e-9, 'RelTol', 1e-7, 'MaxStep', 0.1);
    
    t_traj = linspace(0, T, Tsample);
    x_traj = zeros(Xsample, Tsample);
    x0_list = linspace(xrange(1), xrange(2), Xsample);
    
    for i = 1:Xsample
        x0_curr = x0_list(i);
        sol = ode45(f, [0, T], x0_curr, options);
        x_traj(i, :) = deval(sol, t_traj);
    end
    
    save('occ_data.mat', 't_traj', 'x_traj')
end

if KDE
    Nxspacing = 500;
    x_spacing = linspace(xrange(1), xrange(2), Nxspacing);
    
    pdf_ti = zeros(Nxspacing, Tsample);
    for ki = 2:Tsample
        [pdf_ti(:, ki), ~] = ksdensity(x_traj(:, ki), x_spacing,...
            'Support',1e-2*[-1,1]+x_traj([1, Xsample], ki)','BoundaryCorrection','log');
    end
    
    pdf_ti(:, 1) = ones(Nxspacing, 1)/diff(xrange);
    
%     plot(xi, pdf_ti)
end


if PLOT_TRAJ
    figure(1)
    clf
    hold on
    plot(t_traj, x_traj, 'c')
%     for i = 1:Xsample
        
end

if PLOT_OCC    
%     figure(2)
    clf
    hold on 
%     [TT, XX] = meshgrid(t_traj, x_spacing);
%     surf(TT, XX, pdf_ti);

%     t_probe = [1, 50:50:200];
%     t_probe = [1, 20, 100, 200];
    t_probe = [1, 10, 30, 125, 200];

    for ti = 1:length(t_probe)
        i = t_probe(ti);
        plot3(t_traj(i)*ones(Nxspacing,1), x_spacing, pdf_ti(:, i), 'k', 'LineWidth', 3)
    end
    x_stride = 25;
    x_probe = [1, x_stride:x_stride:5000];
    plot(t_traj, x_traj(x_probe, :), 'c')
    
    view(3)
    
    title('Slices of Occupation Measure','fontsize', 18)
    xlabel('time t')
    ylabel('state x')
    zlabel('state slice prob. density')
end
    