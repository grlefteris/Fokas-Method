function [unapp z_side z_corner ang E]=PDtN_Laplace(Npol,Ncol,M,R,rad_z,rot_z,axis_x,axis_y,problem,scale,u)
[z_side z_corner ang  m h]=irpolygon(Npol,Ncol,rot_z,axis_x,axis_y);
A=PLHS_Dirichlet(Npol,Ncol,M,R,m,h,scale);
Am = sum(abs(A),2);
A = A./repmat(Am,1,Npol*Ncol);
if nargin<11
[G lamda]=Prhs3(Npol,M,R,m,h,ang,problem,scale);
else 
%     G=PmRHS4(Npol,Ncol,M,R,m,h,scale,u);
    G=Prhs(Npol,Ncol,M,R,m,h,reshape(u,Npol*Ncol,1),scale);
end
G=G./Am;
c=A\G; 
c=real(c);
phi=leg_basis(Ncol);
P=kron(eye(Npol),phi);
unapp=P*c;
unapp=reshape(unapp,Ncol,Npol);

%===========Display D2N Errors on every polygon's side=====================
if strcmp(problem,'dirichlet') 
for j=1:Npol
    x=linspace(real(z_corner(j)),real(z_corner(j+1)),Ncol);
    y=linspace(imag(z_corner(j)),imag(z_corner(j+1)),Ncol);
    angx=ang(j)-(pi/2);
    q=@(x,y) cos(angx)*(exp(1+(x/scale)).*cos(2+(y/scale))) +...
        sin(angx)*(-exp(1+(x/scale)).*sin(2+(y/scale)));
%     q=@(x,y) cos(angx)*(exp(1+x).*cos(2+y)) +...
%         sin(angx)*(-exp(1+x).*sin(2+y));
    err=max(abs(unapp(:,j)-q(x,y)'));
    fprintf('Error of side %i = %e\n',j,err);
    err=norm(unapp(:,j)-(q(x,y))',inf)/norm(q(x,y),inf);
    fprintf('NormError of side %i = %e\n',j,err);
    E=max(err);
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