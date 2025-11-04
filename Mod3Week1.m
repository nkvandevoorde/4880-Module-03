%% MODULE 03 WEEK 1 &&

%% 3D Ray Propogation %%
clear; close all; clc;

%Conditions
Fl = 6;
d2 = Fl;
%(x, y, z (d1), r, g, b)
ColorSpace = [
               -2, -1, 50, 1, 0, 0;
               0.5, 1, 40, 0, 1, 0;
               1, -2, 30, 0, 0, 1;
               3, -2, 60, 1, 1, 0;
               ];
rays_perPoint = 50; 
thetas_perPoint = deg2rad(20);
R = 1.5; %Aperature lens radius

%Sensor parameters
x_pixels = 200;
y_pixels = 200; 
s_width = 6;
s_height = s_width;
x_sensorMin = -s_width/2;
y_sensorMin = -s_height/2;
x_sensorMax = s_width/2;
y_sensorMax = s_height/2;
x_sensor = linspace(x_sensorMin, x_sensorMax, x_pixels);
y_sensor = linspace(y_sensorMin, y_sensorMax, y_pixels);
[X, Y] = meshgrid(x_sensor, y_sensor);

plot_rays = true;
max_rays = 500; 

%2x2 matrix
M_Lens = [1 0; -(1/Fl) 1];

%Sensor data
sensor_hits = [];   %(x_img, y_img)
sensor_colors = []; %(r, g, b)

plot_count = 0; 

%Ray Propogation and Tracing 
%rng(0) %turn on for repeatability

for p = 1:size(ColorSpace,1)
    
    x_ColorSpace = ColorSpace(p,1);
    y_ColorSpace = ColorSpace(p,2);
    z_ColorSpace = ColorSpace(p,3);
    ColorSpaceColumns = ColorSpace(p,4:6)'; %RGB columns

    N = rays_perPoint;

    %Generate random rays of propogration to simulate diffuse reflection
    u = rand(1,N);
    v = rand(1,N);
    %angles
    polar = acos(sqrt(u)) * (thetas_perPoint/(pi/2));
    phi = 2*pi*v; %random azimuth angle -> cone orientation
    thx = tan(polar) .* cos(phi); %overall angle of x propogation
    thy = tan(polar) .* sin(phi); %overall angle of y propagtion

    %initial ray positions
    x0 = repmat(x_ColorSpace, 1, N);
    y0 = repmat(y_ColorSpace, 1, N);
    z0 = repmat(z_ColorSpace, 1, N); %distance from lens (cm)

    %RGB per ray
    RGB_ray = repmat(ColorSpaceColumns,1,N);

    %Object -> Lens
    x_lens = x0 + z0 .* thx; %where x propogation intersects lens
    y_lens = y0 + z0 .* thy; %where y propogation intersects lens
    thx_lens = thx; %angle at which x propogation hits lens
    thy_lens = thy; %angle at which y propogation hits lens

    %Aperature check -> if ray is in region then let it through to sensor
    app_region = (x_lens.^2 + y_lens.^2) <= R^2;
    if ~any(app_region)
        continue
    end
    %Only take into account rays through aperature region
    x_lens = x_lens(app_region);
    y_lens = y_lens(app_region);
    thx_lens = thx_lens(app_region);
    thy_lens = thy_lens(app_region);
    ap_columns = RGB_ray(:, app_region);
    ap_N = size(x_lens,2);

    %Ray propogation through lens
    x_atLens = M_Lens * [x_lens; thx_lens]; %[position, angle]
    y_atLens = M_Lens * [y_lens; thy_lens]; %[position, angle]

    x_afterLens = x_atLens(1,:); %x position after lens
    thx_afterLens = x_atLens(2,:); %x angle after lens
    y_afterLens = y_atLens(1,:); %y position after lens
    thy_afterLens = y_atLens(2,:); %y angle after lens

    %Lens -> Img
    x_img = x_afterLens + d2 .* thx_afterLens; %follow x from lens to img
    y_img = y_afterLens + d2 .* thy_afterLens; %follow y from lens to img

    %Img Sensor Coordinated
    sensor_hits = [sensor_hits; x_img(:), y_img(:)]; %fill matrix w/ final x,y
    sensor_colors = [sensor_colors; ap_columns']; %fill matrix w/ final rgb

    %Plot rays
    if plot_rays && plot_count < max_rays
        nplot = min(numel(x_img), max_rays - plot_count);

        figure(1);
        hold on;
        grid on;
        axis equal;
        xlabel('z (cm)');
        ylabel('x (cm)');
        zlabel('y (cm)');
        title('3D Ray Propogation');

        for c = 1:nplot
            %obj -> lens (x)
            plot3([-z_ColorSpace, 0], [x_ColorSpace, x_lens(c)], ...
                [y_ColorSpace, y_lens(c)], 'b--');
            %lens -> Sensor
            plot3([0, d2], [x_afterLens(c), x_img(c)], [y_afterLens(c), ...
                y_img(c)], 'r');
        end
        plot_count = plot_count + nplot; %update plot count
    end
end

%Results
figure; 
scatter(sensor_hits(:,1), sensor_hits(:,2), 6, sensor_colors, 'filled');
xlabel('x_{sensor} (cm)');
ylabel('y_{sensor} (cm)');
title = ('Sensor Positions of Rays');
axis equal;
grid on;

%Aperature region
hold on;
circle = linspace(0, 2*pi, 200);
x_ap = R * cos(circle);
y_ap = R * sin(circle);
plot(x_ap, y_ap, 'w--', 'LineWidth', 1.2); 
legend('Ray hits', 'Aperture');


disp(['Total rays hitting lens: ', num2str(size(sensor_hits,1))]);
