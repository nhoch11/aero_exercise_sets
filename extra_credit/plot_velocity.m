% define filename
text_file = 'hoch_EC.txt';
%text_file = 'cannonball.txt';

% Extract coordinates
data = importdata(text_file, ' ', 1);

time =  data.data(:,1);
u =  data.data(:,2);
v =  data.data(:,3);
w =  data.data(:,4);

x =  data.data(:,8);
y =  data.data(:,9);
z = -data.data(:,10);

% create 3d scatter plot
scatter3(x,y,z);

hold on;

% edit plot labels
title('F-16 Trajectory');
xlabel('X-axis [ft]');
ylabel('Y-axis [ft]');
zlabel('Z-axis [ft]');

% add explosion
lastX = x(end);
lastY = y(end);
lastZ = z(end);



hold off
%aspectRatio = [20,1,1];
%daspect(aspectRatio);
% Show plot
grid on;

figure;
hold on;

% create 3D scatter plot with specified marker face colors
scatter(time, u, 'b', 'filled', 'MarkerFaceColor', 'b');
scatter(time, v, 'r', 'filled', 'MarkerFaceColor', 'r');
scatter(time, w, 'g', 'filled', 'MarkerFaceColor', 'g');

% edit plot labels
title('F-16 Velocity');
xlabel('time [s]');
ylabel('velocity [ft/s]');

legend('u', 'v', 'w', 'location', 'east')

hold off;

% Show plot
grid on;