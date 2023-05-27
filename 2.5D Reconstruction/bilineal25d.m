function [ICoronal, ISagital, IAxial]  = bilineal25d (CT, res, N_ , M_, dx, dy, dz, centro)

% AXES X-Y -> Axial Slice
% X-Z AXES -> Sagittal Slice
% Y-Z AXES -> Coronal Slice

x_c = centro(1);
y_c = centro(2);
z_c = centro(3);

%We choose cut in the center of CT and use reshape to get it in 2D
%ICoronal
img_org1 = reshape(CT(x_c,:,:),size(CT,2),size(CT,3));

%ISagital
img_org2 = reshape(CT(:,y_c,:),size(CT,1),size(CT,3));

%IAxial
img_org3 = reshape(CT(:,:,z_c),size(CT,1),size(CT,2));

%Defino central point
ic1 = y_c;
jc1 = z_c;

ic2 = x_c;
jc2 = z_c;

ic3 = x_c;
jc3 = y_c;

%Multiply the center by its respective differential to obtain it in [mm]
xc1 = ic1*dy;
yc1 = jc1*dz;

xc2 = ic2*dx;
yc2 = jc2*dz;

xc3 = ic3*dx;
yc3 = jc3*dy;

%Position X0 and Y0[mm]
x0_1 = xc1 - (M_/2);
y0_1 = yc1 - (N_/2);

x0_2 = xc2 - (M_/2);
y0_2 = yc2 - (N_/2);

x0_3 = xc3 - (M_/2);
y0_3 = yc3 - (N_/2);

%Define coordinate grid of the image. I create 
%(x,y) pairs for each point on the image that I want to scale.
%I build a single meshgrid because the output dimensions are
%those required by the user regardless of the slice that is made.
[i_, j_] = meshgrid(1:M_/res, 1:N_/res);

%I' y J' with the respective index for each slice
ip1 = (x0_1/dy) + (res/dy)*i_;
jp1 = (y0_1/dz) + (res/dz)*j_;

ip2 = (x0_2/dx) + (res/dx)*i_;
jp2 = (y0_2/dz) + (res/dz)*j_;

ip3 = (x0_3/dx) + (res/dx)*i_;
jp3 = (y0_3/dy) + (res/dy)*j_;

%Calculate floors and ceilings for the bilinear
piso_i1 = floor(ip1);
piso_j1 = floor(jp1);

piso_i2 = floor(ip2);
piso_j2 = floor(jp2);

piso_i3 = floor(ip3);
piso_j3 = floor(jp3);

% Boundary Conditions -> if you go over the size of the image, it will
% impose black pixels

%Coronal images boundary conditions
index1 = piso_i1 < 1;
index2 = piso_i1 > size(img_org1,1) - 1;
index3 = piso_j1 < 1;
index4 = piso_j1 > size(img_org1,2) - 1;

index = ((index1| index2)|(index3|index4));
piso_i1((index)) = 1;
piso_j1((index)) = 1;
ip1(index) = 1;
jp1(index) = 1;

%Sagital images boundary conditions
index12 = piso_i2 < 1;
index22 = piso_i2 > size(img_org2,1) - 1;
index32 = piso_j2 < 1;
index42 = piso_j2 > size(img_org2,2) - 1;

index2 = ((index12| index22)|(index32|index42));
piso_i2(index2) = 1;
piso_j2(index2) = 1;
ip2(index2) = 1;
jp2(index2) = 1;

%Axial images boundary conditions
index13 = piso_i3 < 1;
index23 = piso_i3 > size(img_org3,1) - 1;
index33 = piso_j3 < 1;
index43 = piso_j3 > size(img_org3,2) - 1;

index3 = ((index13| index23)|(index33|index43));
piso_i3(index3) = 1;
piso_j3(index3) = 1;
ip3(index3) = 1;
jp3(index3) = 1;

% a y b for the bilineal
a1 = ip1 - piso_i1;
b1 = jp1 - piso_j1;

a2 = ip2 - piso_i2;
b2 = jp2 - piso_j2;

a3 = ip3 - piso_i3;
b3 = jp3 - piso_j3;

%sub2ind returns the linear indexes equivalent to the row and columnn
%indexes of floor_i and floor_j arrays for the input matrix of size [M,N] 
%This is how I get the indixes of each point I want to access in the for loop 
%of the interpolation instead of having it in matrix form

%Coronal image indexes
%I(i,j)
i_j1 = sub2ind([size(img_org1,1), size(img_org1,2)],piso_i1     , piso_j1);
%I(i+1,j)
i1_j1 = sub2ind([size(img_org1,1), size(img_org1,2)],piso_i1 + 1, piso_j1);
%I(i,j+1)
i_j11 = sub2ind([size(img_org1,1), size(img_org1,2)],piso_i1    , piso_j1 + 1);
%I(i+1,j+1)
i1_j11 = sub2ind([size(img_org1,1), size(img_org1,2)],piso_i1 + 1, piso_j1 + 1);


%Sagital image indexes
%I(i,j)
i_j2 = sub2ind([size(img_org2,1), size(img_org2,2)], piso_i2     , piso_j2);
%I(i+1,j)
i1_j2 = sub2ind([size(img_org2,1), size(img_org2,2)], piso_i2 + 1,piso_j2);
%I(i,j+1)
i_j12 = sub2ind([size(img_org2,1), size(img_org2,2)], piso_i2    , piso_j2 + 1);
%I(i+1,j+1)
i1_j12 = sub2ind([size(img_org2,1), size(img_org2,2)],piso_i2 + 1, piso_j2 + 1);


%Axial image indexes
%I(i,j)
i_j3 = sub2ind([size(img_org3,1), size(img_org3,2)], piso_i3     , piso_j3);
%I(i+1,j)
i1_j3 = sub2ind([size(img_org3,1), size(img_org3,2)], piso_i3 + 1, piso_j3);
%I(i,j+1)
i_j13 = sub2ind([size(img_org3,1), size(img_org3,2)], piso_i3    , piso_j3 + 1);
%I(i+1,j+1)
i1_j13 = sub2ind([size(img_org3,1), size(img_org3,2)], piso_i3 + 1, piso_j3 +1);

%%%%%%%%%%%%%%%%%%%%%%%
%Bilinear interpolation 
%%%%%%%%%%%%%%%%%%%%%%%
ICoronal = ( 1-a1 ).*( 1-b1 ).*img_org1(i_j1)  + ...
           (  a1  ).*( 1-b1 ).*img_org1(i_j11) + ...
           ( 1-a1 ).*(  b1  ).*img_org1(i1_j1) + ...
           (  a1  ).*(  b1  ).*img_org1(i1_j11);
figure;
imshow(imrotate(ICoronal,180));title('Coronal');

ISagital = ( 1-a2 ).*( 1-b2 ).*img_org2(i_j2)  + ...
           (  a2  ).*( 1-b2 ).*img_org2(i_j12) + ...
           ( 1-a2 ).*(  b2  ).*img_org2(i1_j2) + ...
           (  a2  ).*(  b2  ).*img_org2(i1_j12);
figure;
imshow(imrotate(ISagital,180));title('Sagital');

IAxial =   ( 1-a3 ).*( 1-b3 ).*img_org3(i_j3)  + ...
           (  a3  ).*( 1-b3 ).*img_org3(i_j13) + ...
           ( 1-a3 ).*(  b3  ).*img_org3(i1_j3) + ...
           (  a3  ).*(  b3  ).*img_org3(i1_j13);
figure;
imshow(imrotate(IAxial,-90));title('Axial');

end