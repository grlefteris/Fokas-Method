function [A unapp z_side z_corner ang]=PG_D2N_Laplace(Npol,Ncol,M,R,rad_z,rot_z,problem,scale,u)
[z_side z_corner m h ang]=polygon(Npol,Ncol,rad_z,rot_z);
A=PLHS_Dirichlet(Npol,Ncol,M,R,m,h,scale);
Am = sum(abs(A),2);
A = A./repmat(Am,1,Npol*Ncol);
if nargin<9
[G lamda]=Prhs3(Npol,M,R,m,h,ang,problem,scale);
else 
    G=mRHS4(Npol,Ncol,M,R,m,h,u);
end
G=G./Am;
c=A\G; 
c=real(c);
c=reshape(c,Ncol,Npol);
phi=leg_basis(Ncol);
unapp=zeros(Ncol,Npol);
for j=1:Npol
    unapp(:,j)=phi*c(:,j);
end     

%===========Display D2N Errors on every polygon's side=====================
if strcmp(problem,'dirichlet') 
for j=1:Npol
    x=linspace(real(z_corner(j)),real(z_corner(j+1)),Ncol);
    y=linspace(imag(z_corner(j)),imag(z_corner(j+1)),Ncol);
    angx=ang(j)-(pi/2);
    angy=angx+(pi/2);
    q=@(x,y) cos(angx)*(exp(1+(x/scale)).*cos(2+(y/scale))) +...
        sin(angx)*(-exp(1+(x/scale)).*sin(2+(y/scale)));
    err=max(abs(unapp(:,j)-q(x,y)'));
    fprintf('Error of side %i = %e\n',j,err);
    err=norm(unapp(:,j)-(q(x,y))',inf)/norm(q(x,y),inf);
    fprintf('NormError of side %i = %e\n',j,err);
end
else if strcmp(problem,'neumann')
for j=1:Npol
    x=linspace(real(z_corner(j)),real(z_corner(j+1)),Ncol);
    y=linspace(imag(z_corner(j)),imag(z_corner(j+1)),Ncol);
    
    q=@(x,y) exp(1+x).*cos(2+y);
    err=max(abs(unapp(:,j)-q(x,y)'));
    fprintf('Error of side %i = %e\n',j,err);
    err=norm(unapp(:,j)-q(x,y)',inf)/norm(q(x,y),inf);
    fprintf('NormError of side %i = %e\n',j,err);
end
else
        disp('Error: Choose Dirichlet or Neumann Boundary Conditions');
    end
end
disp('\n');
%==========================================================================