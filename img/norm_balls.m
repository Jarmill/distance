pts_inf = [-1 -1 1 1 -1;
           1 -1 -1 1 1];
       
pts_1 = [1 0 -1 0 1;
        0 1 0 -1 0];
    
figure(4)
clf
    
tiledlayout(3, 1)
ax1 = nexttile;
plot(pts_inf(1, :), pts_inf(2, :))
axis square
axis off

ax2 = nexttile;
plot(pts_1(1, :), pts_1(2, :))
axis square
axis off

ax3 = nexttile;
p =3;

th = linspace(0, 2*pi, 200);
xs = cos(th); 
ys = sin(th);

p_norms = (abs(xs).^p + abs(ys).^p).^(1/p);
xp = xs./p_norms;
yp = ys./p_norms;

plot(xp, yp)
axis square
axis off

xlim(1.25*[-1,1])
ylim(1.25*[-1,1])

linkaxes([ax1, ax2, ax3])