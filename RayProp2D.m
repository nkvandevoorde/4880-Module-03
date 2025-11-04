%% MODULE 03 WEEK 1 &&
%% 2D RAY PROPAGATION %%

%Conditions
Fl = 5;
d1 = 30;
d2 = Fl;

%Matrices
M_Fl = [1 d1; 0 1];             %Free space to lens
M_Lens = [1 0; -(1/Fl) 1];      %Lens
M_Fi = [1 d2; 0 1];             %Free space to Image
M = M_Fi * M_Lens * M_Fl;

%2D propagation specs
theta = deg2rad(5);
x = linspace(-2, 2, 5);
y = linspace(-2, 2, 5);
%difuse reflection
rays_perPoint = 10;
thetas_perPoint = deg2rad(5)
R = 10; %appeture radius
sensor_coords = []; %positions on sensor

%Create Visual
figure
hold on;
grid on;
x_label = ('d (cm)');
y_label = ('y (cm)');
title = ('2D Ray Propogation');
axis equal;

%Ray Propagation
rays = numel(y);

for ix = 1:numel(x)
    for iy = 1:numel(y)
        x0_rays = repmat(x(ix), 1, rays_perPoint);
        y0_rays = repmat(y(iy), 1, rays_perPoint);
        thx_rays = theta + (rand(1, rays_perPoint)-0.5) * 2 * thetas_perPoint;
        thy_rays = theta + (rand(1, rays_perPoint)-0.5) * 2 * thetas_perPoint;
        
        for k = 1:rays_perPoint
            x0 = x0_rays(k);
            y0 = y0_rays(k);
            thX0 = thx_rays(k);
            thY0 = thy_rays(k);

            %Free space to lens
            x_rays_lens = M_Fl * [x0; thX0]; %Where x rays hit the lens
            y_rays_lens = M_Fl * [y0; thY0]; % Where y rays hit the lens

            %Aperature (2D)
            if sqrt(x_rays_lens(1)^2 + y_rays_lens(1)^2) > R
                continue;
            end

        %Lens shift
        x_rays_after_lens = M_Lens * x_rays_lens; %How lens shifts x ray
        y_rays_after_lens = M_Lens * y_rays_lens; %How lens shifts y ray

        %Lens to image
        x_ray_img = M_Fi * x_rays_after_lens;
        y_ray_img = M_Fi * y_rays_after_lens;
    
        %Plot rays
        plot([0, d1], [x0, x_rays_lens(1)], 'b--');
        plot([d1, d1 + d2], [x_rays_after_lens(1), x_ray_img(1)], 'b');
        plot([0, d1], [y0, y_rays_lens(1)], 'g--');
        plot([d1, d1 + d2], [y_rays_after_lens(1), x_ray_img(1)], 'g');

        %Find convergence (2D)
        sensor_coords(end+1,:) = [x_ray_img(1), y_ray_img(1)];
        end
    end
end

%Results 
%lens and img position
xline(d1, 'k--', 'LineWidth', 1, 'Label', 'Lens');
xline(d1 + d2, 'g--', 'LineWidth', 1, 'Label', 'Img');
legend('Before Lens', 'After Lens', 'Location', 'best');
hold off;

figure;
scatter(sensor_coords(:,1), sensor_coords(:,2), 30, 'filled');
xlabel('x_{sensor} (cm)');
ylabel('y_{sensor} (cm)');
title = ('Sensor Positions of Rays');
grid on;

%Propagation Matrix
disp('Matrix Multiplication Obj -> Img: ');
disp(M);
