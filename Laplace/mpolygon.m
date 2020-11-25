function [z_side z_corner m h ang ]=mpolygon(Npol,Ncol,rad_z,rot_z)
%==================Polygon's Details=======================================
% t=linspace(-1,1,Ncol);                      %parameter of the side of the polygon
% tn=abs(t(1)-t(2));
% clear t;
t=zeros(1,Ncol);
tn=2/(Ncol-1);
t((Ncol+1)/2)=0;

%==========================================================================

%==========Corners of the Polygon==========================================
z_corner=zeros(Npol+1,1);
for j=1:Npol
    z_corner(j)=rad_z*(cos((2*(j-1)*(pi/Npol))-rot_z)+...
        1i*sin((2*(j-1)*(pi/Npol))-rot_z));
end
    z_corner(Npol+1)=z_corner(1);
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
        z_side(i,j)=m(j)+t(i)*h(j);
    end
    %==============================    
    
end
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