function [z_side z_corner ang m h]=irpolygon(Npol,Ncol,rot_z,axis_x,axis_y)
%==================Polygon's Details=======================================
s=linspace(-1,1,Ncol);
%==========================================================================
%==========================================================================
%=====================Iregular Polygon's Corners===========================
z_corner=zeros(Npol+1,1);

for j=1:Npol
    z_corner(j)=complex(axis_x(j),axis_y(j));
end
z_corner(Npol+1)=z_corner(1);
z_corner=exp(1i*(rot_z))*z_corner;
%==========================================================================                             
%==========Polygon's Angles================================================
ang=zeros(Npol,1);
for j=1:Npol-1
    ang(j)=atan2(imag(z_corner(j+1)-z_corner(j)),...
        real(z_corner(j+1)-z_corner(j)));
end
    ang(Npol)=atan2(imag(z_corner(1)-z_corner(Npol)),...
        real(z_corner(1)-z_corner(Npol)));
%==========================================================================    
%=========Parameterization of the Polygon==================================
z_side=zeros(Ncol,Npol);
m=zeros(Npol,1);
h=zeros(Npol,1);
for j=1:Npol
    
    m(j)=(z_corner(j+1)+z_corner(j))/2;
    h(j)=(z_corner(j+1)-z_corner(j))/2;
    
       
    %parameterize every side=======
    for i=1:Ncol
        z_side(i,j)=m(j)+s(i)*h(j);
    end
    %==============================    
    
end
%==========================================================================
