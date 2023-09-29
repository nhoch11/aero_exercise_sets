% define filename
text_file = 'cannonball.txt';

% Extract coordinates
data = importdata(text_file, ' ', 1);

x =  data.data(:,8);
y =  data.data(:,9);
z = -data.data(:,10);

% create 3d scatter plot
scatter3(x,y,z);

hold on;

% edit plot labels
title('Cannonball Trajectory');
xlabel('X-axis [ft]');
ylabel('Y-axis [ft]');
zlabel('Z-axis [ft]');

% add explosion
lastX = x(end);
lastY = y(end);
lastZ = z(end);


explosionImage = imread('explosion.png');
cannonInput = imread('C:\USU 2023 b Fall\MAE 6510 Aircraft Simulation\aero_exercise_sets\Exercise Set 4\4.3_spin\cannon.png');

cannonImage = imresize(cannonInput, [size(cannonInput, 1), size(cannonInput, 1)*20]);

imsurf(explosionImage, [lastX, 7+lastY, 5], [-1.0,0.0,0.0], [0.0,-1.0,0.0], 0.02)
imsurf(cannonImage, [-300, 0, z(1)+5], [0.0,-1.0,0.0], [0.5,0.0,0.0], 0.02)

hold off
aspectRatio = [20,1,1];
daspect(aspectRatio);
% Show plot
grid on;
