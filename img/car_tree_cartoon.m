figure(808)
clf

%colors
brown = [100 30 35]./255;
black = [0 0 0];
green = [50 150 70]./255;
purple = [107 76 154]./255;
light_purple = [189, 171, 217]/255;


%% draw the tree
c_tree = [0,0];

leaf_width = 1;
leaf_height = 1.5;
leaf = [-leaf_width 0 leaf_width; 0 leaf_height 0];

leaf_stagger = leaf_height*[0, 0.5, 1];
for i = 1:length(leaf_stagger)
patch(leaf(1, :)+c_tree(1), leaf(2, :) + leaf_stagger(i)+c_tree(2), green, 'edgecolor', 'none');
end

trunk_width = 0.4;
trunk_height = 0.6;
trunk = [[-1, 1, 1, -1]*trunk_width;[0, 0, -1, -1]*(trunk_height)];

patch(trunk(1, :)+c_tree(1), trunk(2, :)+c_tree(2), brown, 'edgecolor', 'none');




%% draw the car
c_car = [-3; -2];


%chassis
chas_bot = 1.5;
chas_top = 1.5;
chas_h = 0.2;
chas = [-chas_bot chas_bot chas_top -chas_top;
        0 0 chas_h chas_h] + c_car;
    
patch(chas(1, :), chas(2, :), purple, 'edgecolor', 'none')    

cabin_bot = 1.25;
cabin_top = 0.5;
cabin_h = 1;
cabin = [-cabin_bot cabin_bot cabin_top -cabin_top;
        0 0 cabin_h cabin_h] + c_car + [0; chas_h];
patch(cabin(1, :), cabin(2, :), purple, 'edgecolor', 'none')

%wheels
Nth = 300;
theta = linspace(0, 2*pi, Nth);
circ = [cos(theta); sin(theta)];

Rwheel = 0.2;
dist_wheel = 0.7;
wheel = circ*Rwheel + c_car;


patch(wheel(1, :) + dist_wheel, wheel(2, :), black)
patch(wheel(1, :) - dist_wheel, wheel(2, :), black)


pbaspect([diff(xlim), diff(ylim), 1])