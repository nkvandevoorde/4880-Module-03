%% Parallel rays converging to one point %%
clear; close all; clc; 

%System
d1 = 120;
Fl = 6;
%d2 = (Fl*d1)/(d1-Fl);
d2 = Fl;
thetaX = deg2rad(0);

rowsPositionMatrix = [1, thetaX; 
                      0, thetaX];
PositionMatrix = rowsPositionMatrix';

%Propogation Matrix 
Md1 = [1, d1; 0, 1];
M_Lens = [1, 0; -1/Fl, 1];
Md2 = [1, d2; 0, 1];
M = Md2 * M_Lens * Md1;

%Ray's final postion
Ray_Final = M * PositionMatrix
