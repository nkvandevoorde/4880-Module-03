%% MODULE 03 WEEK 2 &&

%% 3D Ray Propogation %%
clear; close all; clc;

d1_all = linspace(40, 60, 5);   %multiple lens positions
numD1 = numel(d1_all);
Fl = 6;
d2 = -Fl; %focuses at infinity
thetaX = deg2rad(5);
thetaY = deg2rad(5);

%Object Matrix
rowsPositionMatrix = [2, thetaX, 2, thetaY, 1, 0, 0; %objects' position in free space
                      1, thetaX, 1, thetaY, 0, 1, 0;
                      0, thetaX, 0, thetaY, 0, 0, 1;
                      -1, thetaX, -1, thetaY, 0, 1, 1;
                      -2, thetaX, -2, thetaY, 1, 1, 0];
PositionMatrix = rowsPositionMatrix';   %each column is an object
numObj = size(PositionMatrix,2);        %Number of objects
colors = PositionMatrix(5:7,:)';        %Color of object



%Propogation Matricies
M_Lens = [1 0; -(1/Fl) 1];
M_Fi   = [1 d2; 0 1];

%Calculate different obj -> lens matrix for each d1
M_Fl_all = zeros(2,2,numD1);
M_Fl_all(1,1,:) = 1;
M_Fl_all(1,2,:) = d1_all;
M_Fl_all(2,2,:) = 1;

%Propogation Matrix
M = pagemtimes(M_Fi * M_Lens, M_Fl_all);

%2D Propogatiaon Matricies
M_2D = zeros(4,4,numD1);
M_2D(1:2,1:2,:) = M; %x - direction 
M_2D(3:4,3:4,:) = M; %y - direction

% Input rays
Positions = PositionMatrix(1:4,:);  %initial positions

%Propagate rays for all lens positions
Obj_to_Img = pagemtimes(M_2D, Positions);   %[4 x numObj x numD1]
x_img = squeeze(Obj_to_Img(1,:,:));         %[numObj x numD1]
y_img = squeeze(Obj_to_Img(3,:,:));
thetaX = squeeze(Obj_to_Img(2,:,:));
thetaY = squeeze(Obj_to_Img(4,:,:));

%Refocus Parameters
z_focus = 50;  %target refocus plane
x_refocus = x_img + thetaX .* (z_focus - d2); %rays' distance from refocus plane
y_refocus = y_img + thetaY .* (z_focus - d2);

%Plot all rays
obj_idx = [1 2 3 4 5]; %objects
colors_plot = colors(obj_idx, :);

%Planes
z_obj = 0;
z_lenses = d1_all;
z_img = d2;
z_focus = 50;

figure; hold on;

for i = 1:length(obj_idx)
    obj = obj_idx(i);

    for j = 1:numD1
        %Propagation matrices for respective d1
        M_obj_to_lens = [1 d1_all(j); 0 1];
        M_lens = [1 0; -(1/Fl) 1];
        M_lens_to_img = [1 d2; 0 1];

        %Initial ray
        x0 = PositionMatrix(1, obj);
        theta0 = PositionMatrix(2, obj);

        %Object -> Lens
        ray_obj = linspace(z_obj, z_lenses(j), 20);
        x_obj = x0 + theta0 * (ray_obj - z_obj);

        %Apply lens
        ray_state = M_lens * (M_obj_to_lens * [x0; theta0]);
        x_lens = ray_state(1);
        theta_lens = ray_state(2);

        %Lens -> Image
        ray_img = linspace(z_lenses(j), z_img, 20);
        x_img_path = x_lens + theta_lens * (ray_img - z_lenses(j));

        %Extend to refocus plane
        ray_ref = linspace(z_img, z_focus, 20);
        x_ref_path = x_img_path(end) + theta_lens * (ray_ref - z_img);

        %Plot all segments
        plot([ray_obj, ray_img, ray_ref], ...
             [x_obj, x_img_path, x_ref_path], ...
             'Color', colors_plot(i,:), 'LineWidth', 1.2);

        %Ray intersection points
        plot(z_lenses(j), x_lens, 'ko', 'MarkerSize', 4);
    end
end

%Mark planes
xlim([min(d1_all)-10, z_focus+10]);
xlabel('z-axis (optical axis)');
ylabel('x position');
title('Ray Propagation through Plenoptic Setup');
grid on;
yline(0,'k--','LineWidth',0.5);
xline(0,'b--','Object Plane');
xline(mean(d1_all),'g--','Lens Plane');
xline(d2,'r--','Image Plane');
xline(z_focus,'m--','Refocus Plane');
legend({'Ray paths','Lens positions'},'Location','best');
hold off;


%Flatten for plotting
x_plot = -x_refocus(:);       %[numObj*numD1 x 1]
y_plot = -y_refocus(:);       

%Repeat colors for each d1 correctly
color_plot = repmat(colors, numD1, 1);  %[numObj*numD1 x 3]

%Plot
figure;
scatter(x_plot, y_plot, 60, color_plot, 'filled');
xlabel('x_{image}');
ylabel('y_{image}');
title(sprintf('Refocused Image at z = %.1f', z_focus));
axis equal;
grid on;

%Combine all lens positions for each object by averaging (refocus)
x_refocus_mean = mean(x_refocus, 2);  % average across d1 dimension
y_refocus_mean = mean(y_refocus, 2);

%Plot refocused for theta > 0
figure;
scatter(x_refocus_mean, y_refocus_mean, 100, colors, 'filled');
xlabel('x_{image}');
ylabel('y_{image}');
title(sprintf('Refocused Image (averaged over %d lens positions)', numD1));
axis equal;
grid on;
