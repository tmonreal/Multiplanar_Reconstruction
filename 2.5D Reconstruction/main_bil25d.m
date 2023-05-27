% Read DICOM images
S=dicomread('Ima00001');
[M, N] = size(S);
CT=zeros(M,N,140); %We have 140 slices

for i=1:140
    filename=sprintf('Ima%05d',i);
    I=dicomread(filename);
    Iu=double(I)/double(max(I(:)));
    CT(:,:,i) = Iu;  
end

%We define parameters to perform the interpolation
%Desired target image resolution
res = 0.5; 
%Target image dimensions [mm]
N_ = 256; 
M_ = 256; 

%Center for the three slices [pixels]
x_c = 256; 
y_c = 256; 
z_c = 70; 
centro = [x_c, y_c, z_c];

%Measurements of the original images [mm/pix]
dx = 0.459;
dy = 0.459; 
dz = 1.2; 

% Bilinear 2.5d. Receives as parameters:
% - CT: 140 axial tomography slices of 512x512 pixels
% - res: desired output image resolution
% - N_ and M_: measurement in mm of the destination image
% - dx, dy, dz: original image differentials in mm/pix
% - center = [x_c, y_c, z_c]: coordinates of the point chosen as the center

[ICoronal, ISagital, IAxial]  = bilineal25d(CT, res, N_ , M_, dx, dy, dz, centro);

% We can see that the vectorized bilinear is much faster than
% doing it with double for loops. However, the image has outlines
% that are quite pixelated.