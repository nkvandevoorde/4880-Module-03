%% MODULE 03 WEEK 2 &&

%% 3D Ray Propogation %%
clear; close all; clc;

d1 = 50;
Fl = 6;
d2 = -Fl; %focuses at infinity
thetaX = deg2rad(0);
thetaY = deg2rad(0);

%Object Matrix
rowsPositionMatrix = [2, thetaX, 2, thetaY, 1, 0, 0; %objects' position in free space
                      1, thetaX, 1, thetaY, 0, 1, 0;
                      0, thetaX, 0, thetaY, 0, 0, 1;
                      -1, thetaX, -1, thetaY, 0, 1, 1;
                      -2, thetaX, -2, thetaY, 1, 1, 0];
PositionMatrix = rowsPositionMatrix'; %each column is an object

%RGB values
r = PositionMatrix(5,:);
g = PositionMatrix(6,:);
b = PositionMatrix(7,:);

%Color Matrix
colors = [r', g', b'];

%Propogation Matricies
M_Fl = [1 d1; 0, 1];
M_Lens = [1 0; -(1/Fl) 1];      %Lens
M_Fi = [1 d2; 0 1];             %Free space to Image
M = M_Fi * M_Lens * M_Fl;
zeroPadding = [0, 0; 0, 0];
M_FL_2D = [M_Fl, zeroPadding; zeroPadding, M_Fl];
M_2D = [M, zeroPadding; zeroPadding, M];

%%Ray Propogation Obj -> Lens -> Lens
ObjectstoLensMatrix = M_FL_2D * PositionMatrix(1:4, :)
ObjectstoImgMatrix = M_2D * PositionMatrix(1:4, :)

%Initial Positions
x0 = PositionMatrix(1,:);
thx0 = PositionMatrix(2,:);
y0 = PositionMatrix(3,:);
thy0 = PositionMatrix(4,:);

%Lens Positions
x_lens = ObjectstoLensMatrix(1,:);
y_lens = ObjectstoLensMatrix(3,:);

%Image Positions
x_img = ObjectstoImgMatrix(1,:);
y_img = ObjectstoImgMatrix(3,:);

%Plots

%Object -> Lens
figure;
scatter(x_lens, y_lens, 50, 'filled')
title('Ray Intersections at Lens Plane');
xlabel('x'); 
ylabel('y'); 
axis equal; 
grid on;

%Object -> Image
figure;
plot(x0, y0, 'ro-', 'DisplayName', 'Objects');
hold on;
plot(x_lens, y_lens, 'bo-', 'DisplayName', 'Lens');
plot(x_img, y_img, 'go-', 'DisplayName', 'Image');
xlabel('X Position');
ylabel('Y Position');
title('Ray Propagation from Objects to Lens and Image');
legend show;
grid on;

%Reconstructed Image
figure;
scatter(x_img, y_img, 50, colors, 'filled');
title('Reconstructed Image Positions');
xlabel('X Position');
ylabel('Y Position');
axis equal;
grid on;

