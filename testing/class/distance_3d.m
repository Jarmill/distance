BOX = 4;
Ns = 20;

axspace = linspace(-BOX, BOX, Ns);
[XX, YY, ZZ] = meshgrid(axspace, axspace, axspace);

xx = reshape(XX, [], 1);
yy = reshape(YY, [], 1); 
zz = reshape(ZZ, [], 1);

f = @(x, y, z) [y - z.^3, z - y.^2, -x - 2*y - z + y.^2]';
opts = odeset('Events', @(t, x) box_event(t, x, BOX));

figure(1)
clf
hold on
Tmax = 4;
for i = 1:4:Ns
    for j = 1:4:Ns
        for k = 1:4:Ns
            x0 = axspace([i, j, k]);
            [curr_t, curr_x] = ode45(@(t, x) f(x(1), x(2), x(3)), [0, Tmax], x0, opts);
            plot3(curr_x(1, :), curr_x(2, :), curr_x(3, :), 'c')
        end
    end
end


% flow = f(xx, yy, zz);
% dX = reshape(flow(:, 1), Ns, Ns, Ns);
% dY = reshape(flow(:, 2), Ns, Ns, Ns);
% dZ = reshape(flow(:, 3), Ns, Ns, Ns);
% XYZ = stream3(XX, YY, ZZ, dX, dY, dZ, XX, YY, ZZ);
% 
% figure(1)
% clf
% streamline(XYZ)
