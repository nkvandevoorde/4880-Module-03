%% Multiple rays per point %%
clear; close all; clc; 

%System
d1 = 800;
Fl = 6;
%d2 = (Fl*d1)/(d1-Fl);
d2 = Fl;
thetaX0 = deg2rad(linspace(0,180,90));
thetaY0 = deg2rad(linspace(0,180,90));
thetaX1 = deg2rad(linspace(0,180,90));
thetaY1 = deg2rad(linspace(0,180,90));
thetaX2 = deg2rad(linspace(0,180,90));
thetaY2 = deg2rad(linspace(0,180,90));

%Holy Matricies
%--
%(x, y) = (0,0) while theta changes -- simulates cone on rays
%Object 0's x and thetaX
point0X = [zeros(1,90);
            thetaX0];
%Object 0's y and thetaY
point0Y = [zeros(1,90);
            thetaY0];
%--

%--
%(x, y) = (1,1) while theta changes -- simulates cone on rays
%Object 1's x and thetaX
point1X = [ones(1,90);
            thetaX1];
%Object 1's y and thetaY
point1Y = [ones(1,90);
            thetaY1];
%--

%--
%(x, y) = (-1,-1) while theta changes -- simulates cone on rays
%Object -1's x and thetaX
pointNeg1X = [-1*ones(1,90);
            thetaX2];

%Object -1's y and thetaY
pointNeg1Y = [-1*ones(1,90);
            thetaY2];
%-- 

%Positions matrix
point0 = [point0X; point0Y]; %all rays for point (0,0)
point1 = [point1X; point1Y]; %all rays for point (1,1)
pointNeg1 = [pointNeg1X; pointNeg1Y]; %all rays for point (-1,-1)
PositionsMatrix = [point0, point1, pointNeg1]; %All rays for all points

%Propogation Matrix 
Md1 = [1, d1; 0, 1]; %obj -> lens
M_Lens = [1, 0; -1/Fl, 1]; %lens shift
Md2 = [1, d2; 0, 1]; %lens -> sensor
M = Md2 * M_Lens * Md1; %propogation matrix
M_padded = [M, zeros(2); zeros(2), M]; %zero padding to satisfy dimensions


%Ray's final postion
Ray_Final = M_padded * PositionsMatrix; %Where each ray hits sensor

%debugging
%Final0 = M_padded * point0; %Where object 0 hits sensor
%Final1 = M_padded * point1; %Where object 1 hits sensor
%FinalNeg1 = M_padded * pointNeg1; %Where object -1 hits sensor

%Binning and Refocusing
%Move d1 
d700 = 700;
d300 = 300;
d250 = 250;

%700 Matrix
Md700 = [1, d700; 0, 1]; %obj -> lens
M_Lens = [1, 0; -1/Fl, 1]; %lens shift
Md2 = [1, d2; 0, 1]; %lens -> sensor
M = Md2 * M_Lens * Md700; %propogation matrix
M_padded700 = [M, zeros(2); zeros(2), M]; %zero padding to satisfy dimensions

%300 Matrix
Md300 = [1, d300; 0, 1]; %obj -> lens
M_Lens = [1, 0; -1/Fl, 1]; %lens shift
Md2 = [1, d2; 0, 1]; %lens -> sensor
M = Md2 * M_Lens * Md300; %propogation matrix
M_padded300 = [M, zeros(2); zeros(2), M]; %zero padding to satisfy dimensions

%250 Matrix
Md250 = [1, d250; 0, 1]; %obj -> lens
M_Lens = [1, 0; -1/Fl, 1]; %lens shift
Md2 = [1, d2; 0, 1]; %lens -> sensor
M = Md2 * M_Lens * Md250; %propogation matrix
M_padded250 = [M, zeros(2); zeros(2), M]; %zero padding to satisfy dimensions

%Calculate the final positions for each distance matrix
Ray_Final_700 = M_padded700 * PositionsMatrix;
Ray_Final_300 = M_padded300 * PositionsMatrix;
Ray_Final_250 = M_padded250 * PositionsMatrix;

%Combine all rays for blurred image
R_blur = [Ray_Final_250, Ray_Final_300, Ray_Final_700];

% Extract x and y sensor positions
x_sensor = R_blur(1,:); % x position on sensor
y_sensor = R_blur(3,:); % y position on sensor

%Define sensor grid
pixel_size = 0.05;
sensor_half_width = 10;
sensor_half_height = 10;
x_min = -sensor_half_width; x_max = sensor_half_width;
y_min = -sensor_half_height; y_max = sensor_half_height;

x_edges = x_min:pixel_size:x_max;
y_edges = y_min:pixel_size:y_max;
x_centers = (x_edges(1:end-1) + x_edges(2:end))/2;
y_centers = (y_edges(1:end-1) + y_edges(2:end))/2;

%Bin rays
ix = floor((x_sensor - x_min)/pixel_size) + 1;
iy = floor((y_sensor - y_min)/pixel_size) + 1;
valid = ix>=1 & ix<=length(x_centers) & iy>=1 & iy<=length(y_centers);
ix = ix(valid); iy = iy(valid);
sensor_counts = accumarray([ix(:) iy(:)],1,[length(x_centers) length(y_centers)]);
sensor_image = sensor_counts';

%Plot blurred image
figure;
imagesc(x_centers, y_centers, sensor_image);
axis image; xlabel('x (mm)'); ylabel('y (mm)');
title('Blurred Sensor Image');
colorbar;

%Manual Refocus for Each Lens Position
d2 = Fl;  %distance from lens to sensor

%Sensor grid for plotting
pixel_size = 0.05;
sensor_half_width = 10;
sensor_half_height = 10;
x_min = -sensor_half_width; x_max = sensor_half_width;
y_min = -sensor_half_height; y_max = sensor_half_height;
x_edges = x_min:pixel_size:x_max;
y_edges = y_min:pixel_size:y_max;
x_centers = (x_edges(1:end-1) + x_edges(2:end))/2;
y_centers = (y_edges(1:end-1) + y_edges(2:end))/2;

%Refocus for d1 = 250
R = Ray_Final_250;
x_sensor = R(1,:); 
theta_x = R(2,:);
y_sensor = R(3,:); 
theta_y = R(4,:);

%Back-project rays to object plane at d1 = 250
x_ref_250 = x_sensor + theta_x * (d2);
y_ref_250 = y_sensor + theta_y * (d2);

%Bin rays to pixels
ix = floor((x_ref_250 - x_min)/pixel_size) + 1;
iy = floor((y_ref_250 - y_min)/pixel_size) + 1;
valid = ix>=1 & ix<=length(x_centers) & iy>=1 & iy<=length(y_centers);
ix = ix(valid); iy = iy(valid);
ref_counts_250 = accumarray([ix(:) iy(:)], 1, [length(x_centers) length(y_centers)]);
ref_image_250 = ref_counts_250';

%Refocus for d1 = 300
R = Ray_Final_300;
x_sensor = R(1,:); 
theta_x = R(2,:);
y_sensor = R(3,:); 
theta_y = R(4,:);

x_ref_300 = x_sensor + theta_x * (d2);
y_ref_300 = y_sensor + theta_y * (d2);

ix = floor((x_ref_300 - x_min)/pixel_size) + 1;
iy = floor((y_ref_300 - y_min)/pixel_size) + 1;
valid = ix>=1 & ix<=length(x_centers) & iy>=1 & iy<=length(y_centers);
ix = ix(valid); iy = iy(valid);
ref_counts_300 = accumarray([ix(:) iy(:)], 1, [length(x_centers) length(y_centers)]);
ref_image_300 = ref_counts_300';

%Refocus for d1 = 700
R = Ray_Final_700;
x_sensor = R(1,:); 
theta_x = R(2,:);
y_sensor = R(3,:); 
theta_y = R(4,:);

x_ref_700 = x_sensor + theta_x * (d2);
y_ref_700 = y_sensor + theta_y * (d2);

ix = floor((x_ref_700 - x_min)/pixel_size) + 1;
iy = floor((y_ref_700 - y_min)/pixel_size) + 1;
valid = ix>=1 & ix<=length(x_centers) & iy>=1 & iy<=length(y_centers);
ix = ix(valid); iy = iy(valid);
ref_counts_700 = accumarray([ix(:) iy(:)], 1, [length(x_centers) length(y_centers)]);
ref_image_700 = ref_counts_700';

%Plot all refocused images
figure;
subplot(1,3,1);
imagesc(x_centers, y_centers, ref_image_250);
axis image; xlabel('x'); ylabel('y'); title('Refocus at d1 = 250');

subplot(1,3,2);
imagesc(x_centers, y_centers, ref_image_300);
axis image; xlabel('x'); ylabel('y'); title('Refocus at d1 = 300');

subplot(1,3,3);
imagesc(x_centers, y_centers, ref_image_700);
axis image; xlabel('x'); ylabel('y'); title('Refocus at d1 = 700');