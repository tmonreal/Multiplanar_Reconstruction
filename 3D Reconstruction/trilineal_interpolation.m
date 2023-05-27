function pixelInterpolado=trilineal_interpolation(X,CT)

%Point to interpolate
xc = X(1);
yc = X(2);
zc = X(3);

%Calculate floors and ceilings for interpolation
piso_i = floor(xc);
piso_j = floor(yc);
piso_k = floor(zc);

%Boundary conditions -> if it falls outside of the volume, we assign 0
    if( (piso_i < 1) || (piso_i > size(CT,1)-1) )
        xc = -1;   
    end
    
    if( (piso_j < 1) || (piso_j > size(CT,2)-1) )
        yc = -1;
    end
    
    if( (piso_k < 1) || (piso_k > size(CT,3)-1) )
        zc = -1;
    end
 
if ( xc == -1 || yc == -1 || zc == -1 )
        pixelInterpolado = 0;
        
else   
% a, b y c para la trilineal
a = xc - piso_i;
b = yc - piso_j;
c = zc - piso_k;


pixelInterpolado =(1-a).*(1-b).*(1-c).* CT(piso_i    , piso_j    , piso_k    ) + ...
                  ( a ).*(1-b).*(1-c).* CT(piso_i + 1, piso_j    , piso_k)     + ...
                  ( a ).*(1-b).*( c ).* CT(piso_i + 1, piso_j    , piso_k + 1) + ...
                  (1-a).*(1-b).*( c ).* CT(piso_i    , piso_j    , piso_k + 1) + ...
                  (1-a).*( b ).*( c ).* CT(piso_i    , piso_j + 1, piso_k + 1) + ...
                  (1-a).*( b ).*(1-c).* CT(piso_i    , piso_j + 1, piso_k    ) + ...
                  ( a ).*( b ).*(1-c).* CT(piso_i + 1, piso_j + 1, piso_k    ) + ...
                  ( a ).*( b ).*( c ).* CT(piso_i + 1, piso_j + 1, piso_k + 1);
end
              
end




