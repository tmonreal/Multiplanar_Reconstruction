%Load images
S=dicomread('Ima00001');
[M, N] = size(S);
CT=zeros(M,N,140); %We have 140 slices

for i=1:140
    filename=sprintf('Ima%05d',i);
    I=dicomread(filename);
    Iu=double(I)/double(max(I(:)));
    CT(:,:,i) = Iu;   
end

%Canonical basis
i=[1 0 0];
j=[0 1 0];
k=[0 0 1];

%Center of pixels
x_c=256; 
y_c=256; 
z_c=70; 

c=[x_c;y_c;z_c];

%Output image dimensions [pixels] = [512, 512]
M_=512;
N_=512;

%You can modify the output resolution by changing the parameter res
%res = [mm/pix]
res=0.5;

%Input image measurements
dx = 0.459; %[mm/pixel]
dy = 0.459; %[mm/pixel]
dz = 1.2;   %[mm/pixel]
dp=[1/dx;1/dy;1/dz]; %[pixel/mm]

u = i';
v = j';
n = k';

%Output size in [mm]
h=M_*res; 
w=N_*res;

%Versor to rotate
versor=v; 

for theta=0:pi/8:2*pi;
matrot=[cos(theta) + versor(1)^2*(1-cos(theta))                     versor(1)*versor(2)*(1-cos(theta))-versor(3)*sin(theta)     versor(1)*versor(3)*(1-cos(theta))+versor(2)*sin(theta);
        versor(2)*versor(1)*(1-cos(theta))+versor(3)*sin(theta)     cos(theta)+(versor(2)^2*(1-cos(theta)))                     versor(2)*versor(3)*(1-cos(theta))-versor(1)*sin(theta);
        versor(3)*versor(1)*(1-cos(theta))-versor(2)*sin(theta)     versor(3)*versor(2)*(1-cos(theta))+versor(1)*sin(theta)     cos(theta)+versor(3)^2*(1-cos(theta))                 ];

%New versors after rotation
u=matrot(:,1);
v=matrot(:,2);
n=matrot(:,3);

%Normalize u
u_=u/norm(u);

%We verify orthogonality and that they are unitary vectors
%The cross product between the versor and another vector (e.g. u and v) returns the perpendicular vector (e.g. n). 
%Then, the cross product between the vector obtained (e.g. n) and one of the original vectors (e.g. u) returns the other original vector (e.g. v)
n_=cross(u_,v);
n_=n_/norm(n_);

v_=cross(n_,u_);
v_=v_/norm(v_);

%We create the point of origin p from the rotated vectors
p = c./dp - (v_).*(h/2)- (n_).*(w/2); %mm

%Resulting image
X=zeros(M_,N_);

%Initialization
i=1; j=1;

for alpha=1:res:h
    for beta=1:res:w
        x = p + alpha*(v_) + beta*(n_); %mm
        x = x.*dp;
        X(i,j)=interp_trilineal(x',CT);
        j=j+1;
    end
j = 1;
i=i+1;
end
    
imshow(imrotate(X,-90));
pause(0.1);

end