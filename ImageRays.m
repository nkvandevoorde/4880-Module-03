%% Image Lightfield Data

%camera specs
Fl = 0.6; %cm
d2 = Fl; 
lens_Height = 3.6; %cm

%image loading
filename = 'FirstImage.tif';
Image1 = imread(filename);
figure
imshow(Image1)
figure
imhist(Image1)

%%Blue-Blue Block
x_coorBB = -1;
y_coorBB =-1;
d1_BB = 75; %cm
halfThetaBB = atan(lens_Height/d1_BB);
thetaBB = linspace(-1*halfTheta0, halfTheta0, 50);

%%Blue-Yellow Block
x_coorBY = 0;
y_coorBY =0;
d1_BY = 131; %cm
halfThetaBY = atan(lens_Height/d1_BY);
thetaBY = linspace(-1*halfTheta0, halfTheta0, 50);

%%Blue-Blue Block
x_coorRG = 1;
y_coorRG =1;
d1_RG = 181; %cm
halfThetaRG = atan(lens_Height/d1_RG);
thetaRG = linspace(-1*halfTheta0, halfTheta0, 50);

%% distance from camera to edge of table is 137.5cm