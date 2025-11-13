%% Multiple rays per point %%
clear; close all; clc; 

%System
d1 = 120;
Fl = 6;
%d2 = (Fl*d1)/(d1-Fl);
d2 = Fl;
thetaX0 = deg2rad(0);
thetaY0 = deg2rad(0);
thetaX1 = deg2rad(45);
thetaY1 = deg2rad(45);
thetaX2 = deg2rad(-45);
thetaY2 = deg2rad(-45);

%Holy Matricies
%For each "object" make a seperate colum for each ray (x and y are the same 
% but theta changes to simulate multiple rays) and start with overall
% postion as two matricies to ensure that x position and y position are in
% seperate columns
%Object 0's x and thetaX
point0X = [0, thetaX1;
          0, thetaX0;
          0, thetaX2;];

%Object 0's y and thetaY
point0Y = [0, thetaY1; 
           0, thetaY0; 
           0, thetaY2;];

%Object 1's x and thetaX
point1X = [1, thetaX1;
          1, thetaX0;
          1, thetaX2;];

%Object 1's y and thetaY
point1Y = [1, thetaY1; 
           1, thetaY0; 
           1, thetaY2;];

%Object -1's x and thetaX
pointNeg1X = [-1, thetaX1;
              -1, thetaX0;
              -1, thetaX2;];

%Object -1's y and thetaY
pointNeg1Y = [-1, thetaY1;
              -1, thetaY0;
              -1, thetaY2;];

point0 = [point0X(1,:)', point0Y(1,:)', point0X(2,:)', point0Y(2,:)', ... 
          point0X(3,:)', point0Y(3,:)']; %Object 0's overall position
point1 = [point1X(1,:)', point1Y(1,:)', point1X(2,:)', point1Y(2,:)', ... 
          point1X(3,:)', point1Y(3,:)']; %Object 1's overall position
pointNeg1 = [pointNeg1X(1,:)', pointNeg1Y(1,:)', pointNeg1X(2,:)', ...
    pointNeg1Y(2,:)', pointNeg1X(3,:)', pointNeg1Y(3,:)']; %Object -1's overall position

PositionsMatrix = [point0, point1, pointNeg1] %All positions

%Propogation Matrix 
Md1 = [1, d1; 0, 1];
M_Lens = [1, 0; -1/Fl, 1];
Md2 = [1, d2; 0, 1];
M = Md2 * M_Lens * Md1;
M_padded = [0 6 0 0;
           -1/6 -19 0 0];


%Ray's final postion
Ray_Final = M * PositionsMatrix %Where each ray hits sensor
Final0 = M * point0; %Where object 0 hits sensor
Final1 = M * point1; %Where object 1 hits sensor
FinalNeg1 = M * pointNeg1; %Where object -1 hits sensor
