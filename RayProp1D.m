%% MODULE 03 WEEK 1 &&
%% 1D RAY PROPAGATION %%

%Conditions
Fl = 5;
d1 = 30;
d2 = Fl; %simulating a camera at d1 = infinity 
theta = deg2rad(0);
y = linspace(-2, 2, 5);

%Matrices
M_Fl = [1 d1; 0 1];             %Free space to lens
M_Lens = [1 0; -(1/Fl) 1];      %Lens
M_Fi = [1 d2; 0 1];             %Free space to Image
M = M_Fi * M_Lens * M_Fl;

%1D propagation specs
propogation_space = linspace(0, d1+d2, 100); %Space from object to image
lens = d1; %lens location
img_pos = d1 + d2; %Where the image converges

%Diffuse reflection
rays_perPoint = 10;
thetas_perPoint = deg2rad(5);

%Initialize sensor array for reconstruction
pixels = 20; %# of pixels on sensor
s_width = 4; %width covered by rays
pixel_edges = linspace(-s_width/2, s_width/2, pixels + 1);
pixel_energy = zeros(1, pixels);  %energy per pixel

%Create Visual
figure
hold on;
grid on;
x_label = ('d (cm)');
y_label = ('y (cm)');
title = ('Ray Propogation');
axis equal;

%Ray Propagation
rays = numel(y);

for i = 1:rays
    y0_rays = repmat(y(i), 1, rays_perPoint);
    theta0_rays = theta + (rand(1,rays_perPoint)-0.5) * 2 * thetas_perPoint;

    for j = 1:rays_perPoint
        y0 = y0_rays(i);
        theta0 = theta0_rays(i);

        %Free space to lens
        l1 = linspace(0, d1, 50);
        y1 = y0 + tan(theta0) * l1;
        plot(l1, y1, 'b--', 'LineWidth', 1.5);
        rays_lens = M_Fl * [y0;theta0]; %Where rays hit the lens

        %Apperature (1D)
        R = 10; 
            if abs(rays_lens(1)) > R
                continue;
            end
    
        %Lens shift
        rays_after_lens = M_Lens * rays_lens; %How lens shifts ray
        y_after_lens = rays_after_lens(1);    %Lens' shift on y
        theta_after_lens = rays_after_lens(2);%Lens' shift on theta

        %Lens to image
        l2 = linspace(d1, d1 + d2, 50);
        % Where ray travels after lens
        y2 = y_after_lens + (l2 - d1) * tan(theta_after_lens);
        plot(l2, y2, 'r', 'LineWidth', 1.5)

        %Find convergence (1D)
        rays_img = M_Fi * rays_after_lens; %where ray hits the sensor
        y_final(j) = rays_img(1);

        %Find activated pixels
        [~, idx] = min(abs(pixel_edges - rays_img(1)));
        if idx > 1 && idx <= pixels
            pixel_energy(idx) = pixel_energy(idx) + 1;
        end
    end
end

%Results 
%lens and img positions
xline(lens, 'k--', 'LineWidth', 1, 'Label', 'Lens');
xline(img_pos, 'g--', 'LineWidth', 1, 'Label', 'Img');
legend('Before Lens', 'After Lens', 'Location', 'best');
hold off;

%Sensor Reconstruction
figure;
bar(1:pixels, pixel_energy, 'FaceColor', [0.2 0.6 1]);
xlabel('Pixel Index');
ylabel('Ray Energy');
title = ('Sensor Reconstruction');
grid on;

%Propagation Matrix
disp('Matrix Multiplication Obj -> Img: ');
disp(M);
